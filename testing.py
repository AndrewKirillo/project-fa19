import os
import sys
sys.path.append('..')
sys.path.append('../..')
import argparse
import utils
from student_utils import *
from gurobipy import *
import networkx as nx


data_parser("inputs/65_50.in")

num_homes = len(list_of_homes)
num_locations = len(list_of_locations)

# Initializing Graph object
G = nx.Graph()


start_node = list_of_locations.index(starting_car_location)

E = set()

for i in range(len(adjacency_matrix)):
    for j in range(len(adjacency_matrix)):
        if (adjacency_matrix[i][j] != 'x'):
            E.add((i, j))
            G.add_edge(i, j, weight=int(adjacency_matrix[i][j]))

paths, sp_distances = nx.floyd_warshall_predecessor_and_distance(G)



# Initializing capacities
edge_capacity_dict = {}
for edge in E:
    edge_capacity_dict[edge] = 1 

m = Model()


################## VARIABLES ##################

# Edge inclusion variables.
# Edge (i,j) in subtour => edge_vars[i,j] = 1
edge_vars = {}

for edge in E:
    edge_vars[edge] = m.addVar(vtype=GRB.BINARY)


# TA dropoff variables.
# TA t dropped off at location l => ta_vars[t,l] = 1
ta_vars = {}

for t in range(num_homes):
    for l in range(num_locations):
        ta_vars[t,l] = m.addVar(vtype=GRB.BINARY)


# Edge flow variables.
# Edge (i, j) in tour that includes start => flow_vars[i,j] = 1
flow_vars = {}

for edge in E:
    flow_vars[edge] = m.addVar(vtype=GRB.BINARY)
flow_vars[('s', start_node)] = m.addVar(vtype=GRB.INTEGER)
flow_vars[(start_node, 't')] = m.addVar(vtype=GRB.INTEGER)


################# CONSTRAINTS #################

# Max tour length constraint
m.addConstr(quicksum(edge_vars.values()) <= num_locations, "tour_length")

# Cycle enforcement constraint
for v in range(num_locations):
    m.addConstr(quicksum(edge_vars[i,v] if (i, v) in E else 0 for i in range(num_locations)) == quicksum(edge_vars[v,j] if (v, j) in E else 0 for j in range(num_locations)), "cycle_check_node" + str(v))

# One dropoff per TA constraint
for t in range(num_homes):
    m.addConstr(quicksum(ta_vars[t,l] for l in range(num_locations)) == 1, "one_dropoff_ta" + str(t))

# All TAs dropped off
m.addConstr(quicksum(ta_vars.values()) == num_homes, "all_dropped_off")

# TA dropped off on path
for t in range(num_homes):
    for l in range(num_locations):
        m.addConstr((ta_vars[t,l] == 1) >> (max_([edge_vars[l, v] if (l, v) in E else 0 for v in range(num_locations)]) == 1), "ta" + str(t) + "_loc" + str(l) + "_valid_dropoff")


# Flow i,j is less than Capacity i,j
for edge in E:
    m.addConstr(flow_vars[edge] <= edge_capacity_dict[edge], "capacity_not_exceeded")

# Flow passed along path from node to node
for v in range(num_locations)
    m.addConstr(quicksum(flow_vars[i,v] if (i,v) in E else 0 for i in range(num_locations)) == quicksum(flow_vars[v,j] if (v,j) in E else 0 for j in range(num_locations)), "flow_conserved")
m.addConstr(flow_vars['s', start_node] == quicksum(flow_vars[start_node, v] if (start_node, v) in E else 0 for v in range(num_locations)), "flow_start")
m.addConstr(flow_vars[start_node, 't'] == quicksum(flow_vars[u, start_node] if (u, start_node) in E else 0 for u in range(num_locations)), "flow_end")

# Edge included => Flow == 1 == Inclusion Variable

for edge in E:
    m.addConstr(flow_vars[edge] == edge_vars[edge], "inclusion_implies_flow")


################Setting Objective###############
m.setObjective((2.0/3.0)*quicksum(edge_vars[e]*length[e] for e in E) + quicksum(quicksum(ta_vars[t,l]*sp_distances[t][l] for l in range(num_locations)) for t in range(num_homes)), GRB.MINIMIZE)



m.optimize()

