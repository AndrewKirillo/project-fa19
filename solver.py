import os
import sys
sys.path.append('..')
sys.path.append('../..')
import argparse
import utils
from student_utils import *
from gurobipy import *
import networkx as nx

"""
======================================================================
  Complete the following function.
======================================================================
"""

def solve(list_of_locations, list_of_homes, starting_car_location, adjacency_matrix, params=[]):
    """
    Write your algorithm here.
    Input:
        list_of_locations: A list of locations such that node i of the graph corresponds to name at index i of the list
        list_of_homes: A list of homes
        starting_car_location: The name of the starting location for the car
        adjacency_matrix: The adjacency matrix from the input file
    Output:
        A list of locations representing the car path
        A dictionary mapping drop-off location to a list of homes of TAs that got off at that particular location
        NOTE: both outputs should be in terms of indices not the names of the locations themselves
    """

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
        edge_vars[edge] = m.addVar(name="e_"+str(edge), vtype=GRB.BINARY)

    
    # TA dropoff variables.
    # TA t dropped off at location l => ta_vars[t,l] = 1
    ta_vars = {}

    for t in range(num_homes):
        for l in range(num_locations):
            ta_vars[t,l] = m.addVar(name="ta_"+str(t)+"@"+str(list_of_locations[l]), vtype=GRB.BINARY)
    
    
    # Edge flow variables.
    # Edge (i, j) in tour that includes start => flow_vars[i,j] = 1
    flow_vars = {}

    for edge in E:
        flow_vars[edge] = m.addVar(name="flow_"+str(edge), vtype=GRB.BINARY)
    flow_vars[('s', start_node)] = m.addVar(name="flow_source", vtype=GRB.INTEGER)
    flow_vars[(start_node, 't')] = m.addVar(name="flow_target", vtype=GRB.INTEGER)

    
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
            m.addConstr(ta_vars[t,l] <= quicksum([edge_vars[l, v] if (l, v) in E else 0 for v in range(num_locations)]), "ta" + str(t) + "_loc" + str(l) + "_valid_dropoff")


    # Flow i,j is less than Capacity i,j
    for edge in E:
        m.addConstr(flow_vars[edge] <= edge_capacity_dict[edge], "capacity_not_exceeded")

    # Flow passed along path from node to node
    for v in range(num_locations):
        m.addConstr(quicksum(flow_vars[i,v] if (i,v) in E else 0 for i in range(num_locations)) == quicksum(flow_vars[v,j] if (v,j) in E else 0 for j in range(num_locations)), "flow_conserved")
    m.addConstr(flow_vars['s', start_node] == quicksum(flow_vars[start_node, v] if (start_node, v) in E else 0 for v in range(num_locations)), "flow_start")
    m.addConstr(flow_vars[start_node, 't'] == quicksum(flow_vars[u, start_node] if (u, start_node) in E else 0 for u in range(num_locations)), "flow_end")

    # Flow is non-negative
    for edge in E:
        m.addConstr(flow_vars[edge] >= 0, "non-negative flow " + str(edge))
    m.addConstr(flow_vars['s', start_node] >= 0, "non-negative source flow")
    m.addConstr(flow_vars[start_node, 't'] >= 0, "non-negative target flow")
    
    # Edge included => Flow == 1 == Inclusion Variable
    for edge in E:
        m.addConstr(flow_vars[edge] == edge_vars[edge], "inclusion_implies_flow")


    ################Setting Objective###############
    m.setObjective((2.0/3.0)*quicksum(edge_vars[e]*int(adjacency_matrix[e[0]][e[1]]) for e in E) + quicksum(quicksum(ta_vars[t,l]*sp_distances[t][l] for l in range(num_locations)) for t in range(num_homes)), GRB.MINIMIZE)
    
    m.optimize()

    for v in m.getVars():
        vars_to_output(E, adjacency_matrix, start_node, edge_vars)


"""
======================================================================
   No need to change any code below this line
======================================================================
"""

"""
Convert solution with path and dropoff_mapping in terms of indices
and write solution output in terms of names to path_to_file + file_number + '.out'
"""
def convertToFile(path, dropoff_mapping, path_to_file, list_locs):
    string = ''
    for node in path:
        string += list_locs[node] + ' '
    string = string.strip()
    string += '\n'

    dropoffNumber = len(dropoff_mapping.keys())
    string += str(dropoffNumber) + '\n'
    for dropoff in dropoff_mapping.keys():
        strDrop = list_locs[dropoff] + ' '
        for node in dropoff_mapping[dropoff]:
            strDrop += list_locs[node] + ' '
        strDrop = strDrop.strip()
        strDrop += '\n'
        string += strDrop
    utils.write_to_file(path_to_file, string)

def solve_from_file(input_file, output_directory, params=[]):
    print('Processing', input_file)

    input_data = utils.read_file(input_file)
    num_of_locations, num_houses, list_locations, list_houses, starting_car_location, adjacency_matrix = data_parser(input_data)
    car_path, drop_offs = solve(list_locations, list_houses, starting_car_location, adjacency_matrix, params=params)

    basename, filename = os.path.split(input_file)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    output_file = utils.input_to_output(input_file, output_directory)

    convertToFile(car_path, drop_offs, output_file, list_locations)


def solve_all(input_directory, output_directory, params=[]):
    input_files = utils.get_files_with_extension(input_directory, 'in')

    for input_file in input_files:
        solve_from_file(input_file, output_directory, params=params)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Parsing arguments')
    parser.add_argument('--all', action='store_true', help='If specified, the solver is run on all files in the input directory. Else, it is run on just the given input file')
    parser.add_argument('input', type=str, help='The path to the input file or directory')
    parser.add_argument('output_directory', type=str, nargs='?', default='.', help='The path to the directory where the output should be written')
    parser.add_argument('params', nargs=argparse.REMAINDER, help='Extra arguments passed in')
    args = parser.parse_args()
    output_directory = args.output_directory
    if args.all:
        input_directory = args.input
        solve_all(input_directory, output_directory, params=args.params)
    else:
        input_file = args.input
        solve_from_file(input_file, output_directory, params=args.params)
