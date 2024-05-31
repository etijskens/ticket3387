import numpy as np
import pandas as pd
from mesa import Model, Agent
from mesa.time import RandomActivation
from mesa.datacollection import DataCollector

from mesa.batchrunner import batch_run
from multiprocessing import freeze_support

# Load data and other initializations only on the root process
data = pd.read_excel('combined_data_filled_mean.xlsx', engine='openpyxl')
prob_df = pd.read_csv('probability_combinations_subset.csv')
print('len(data) =', len(data))
print('len(prob_df) =', len(prob_df))

class TCellAgent(Agent):
    def __init__(self, unique_id, model, age, gender, status, tcr_affinity, serostatus, initial_state):
        super().__init__(unique_id, model)
        self.age = age
        self.gender = gender
        self.status = status
        self.tcr_affinity = tcr_affinity
        self.serostatus = serostatus
        self.state = initial_state
        self.unique_id = unique_id

    def step(self):
        transition_probabilities = {
            "Naive": ["Effector", "Memory"],
            "Effector": ["Memory", "Dead", "Effector"],
            "Memory": ["Effector", "Dead", "Memory"]
        }
        next_states = transition_probabilities.get(self.state, [])
        if next_states:
            probabilities = []
            for state in next_states:
                try:
                    prob = getattr(self.model, f'p_{self.state.lower()}_{state.lower()}')
                except AttributeError:
                    prob = 0
                probabilities.append(prob)
            total = sum(probabilities)
            if total > 0:
                probabilities = [p / total for p in probabilities]
            else:
                probabilities = [1.0 / len(next_states)] * len(next_states)
            self.state = np.random.choice(next_states, p=probabilities)

class TCellModel(Model):
    def __init__(self, patient_id, p_naive_effector, p_effector_memory, p_memory_effector):
        super().__init__()
        self.data = data[data['Vaccinee'] == patient_id] 
        self.schedule = RandomActivation(self)
        self.p_naive_effector = p_naive_effector
        self.p_naive_memory = 1 - p_naive_effector
        self.p_effector_memory = p_effector_memory
        self.q_effector_dead = 1 - p_effector_memory
        self.p_memory_effector = p_memory_effector
        self.q_memory_dead = 1 - p_memory_effector
        self.datacollector = DataCollector({
            "Naive": lambda m: self.count_state("Naive"),
            "Effector": lambda m: self.count_state("Effector"),
            "Memory": lambda m: self.count_state("Memory"),
            "Dead": lambda m: self.count_state("Dead")
        })
        self.load_data_and_initialize_agents()

    def load_data_and_initialize_agents(self):
        unique_id = 0
        #print('number of rows involved =', len(self.data))
        for _, row in self.data.iterrows():
            age = row['Age']
            gender = row['Gender']
            status = row['Status_2']
            tcr_affinity = row['TCRpredictor']
            serostatus = {
                "CMV": row['CMV_Status'],
                "EBV": row['EBV_Status'],
                "HSV1_2": row['HSV1_2_Status'],
                "HHV6": row['HHV6_Status']
            }
            memory_ratio = row['memory Tconv/ total cd4 t-cells']
            total_cells = 1400
            memory_count = int(total_cells * memory_ratio)
            max_possible_effectors = total_cells - memory_count
            effector_multiplier = 2
            effector_count = min(effector_multiplier * memory_count, max_possible_effectors)
            naive_count = total_cells - memory_count - effector_count
            naive_count = 2
            memory_count = 3
            effector_count = 4
            #print('naive = %s, memory = %s, effector = %s' % (naive_count, memory_count, effector_count))
            for state, count in [("Naive", naive_count), ("Memory", memory_count), ("Effector", effector_count)]:
                for _ in range(count):
                    agent = TCellAgent(unique_id, self, age, gender, status, tcr_affinity, serostatus, state)
                    self.schedule.add(agent)
                    unique_id += 1
        self.num_agents = unique_id
        print('TCellModel initialized, using %s agents' % self.num_agents, flush=True)

    def count_state(self, state):
        return sum(1 for agent in self.schedule.agents if agent.state == state)

    def step(self):
        self.datacollector.collect(self)
        self.schedule.step()

def run_simulation(patient_data, prob_values, steps=365):
    model = TCellModel(patient_data, *prob_values)
    for _ in range(steps):
        model.step()
    results = model.datacollector.get_model_vars_dataframe()
    return results

# Function to prepare the parameters for batch run
def prepare_batch_params(patient_id, prob_df):
    params = []
    for _, row in prob_df.iterrows():
        prob_values = row[['p_naive_to_effector', 'p_effector_to_memory', 'p_memory_to_effector']].values
        param_set = {
            "patient_id": patient_id,
            "p_naive_effector": prob_values[0],
            "p_effector_memory": prob_values[1],
            "p_memory_effector": prob_values[2],
            #"prob_combination": prob_values  # Keep track of probability combination
        }
        params.append(param_set)
    return params

patient_ids = data['Vaccinee'].unique()

if __name__ == "__main__":
    for patient_id in patient_ids:
        #patient_data = data[data['Vaccinee'] == patient_id]
        batch_params = prepare_batch_params(patient_id, prob_df)
        print('len(batch_params) =', len(batch_params))

        freeze_support()
        for params in batch_params:
            print()
            print('params =', params)
            results = batch_run(
                TCellModel,
                parameters=params, #batch_params,
                iterations=1,
                max_steps=3, #365,
                number_processes=None,
                data_collection_period=-1,
                display_progress=True
            )

        # Save results
        for i, (result, param) in enumerate(zip(results, batch_params)):
            prob_values_str = "_".join(map(str, param["prob_combination"]))
            df = pd.DataFrame(result)
            df.to_csv(f"simulation_results_patient{patient_id}_{prob_values_str}.csv")

