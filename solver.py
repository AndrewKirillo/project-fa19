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
    # pass
    home_list = convert_locations_to_indices(list_of_homes, list_of_locations)
    location_list = convert_locations_to_indices(list_of_locations, list_of_locations)
    
    G, message = adjacency_matrix_to_graph(adjacency_matrix)

    #Dropoff mappings dictionary
    dropoff_mappings = {}
#     for home in home_list:
#         dropoff_mapping[home] = [home]

    # Create a minimum spanning_tree of G
    #T = nx.minimum_spanning_tree(G)

    #  Predecessors, distances of the minimum spanning tree
    predecessors, distances = nx.floyd_warshall_predecessor_and_distance(G)

    #list_of_edges = sorted(T.edges(data=False))
    starting_location = list_of_locations.index(starting_car_location)
    curr_location = starting_location
    car_path = [curr_location]
    visited_homes = [] ##aka car_path
    closest_home_path = []
    # visited_homes.append(curr_location)
    not_yet_visited_homes = home_list[:]

    while len(not_yet_visited_homes) != 0:
        closest_home, closest_home_path = find_closest_home_to_current_location(predecessors, distances, visited_homes, not_yet_visited_homes, curr_location)
        visited_homes.append(closest_home)
        not_yet_visited_homes.remove(closest_home)
        curr_location = closest_home
        car_path+=closest_home_path

    # Return home
    path_home = nx.reconstruct_path(curr_location, starting_location, predecessors)
    car_path+=path_home
    
    clean_path = []

    clean_path.append(car_path[0])

    for i in range(1, len(car_path)):
        if car_path[i] != clean_path[-1]:
            clean_path.append(car_path[i])
            
    
    car_path = clean_path[:]
    node_list = car_path[:]
    
    ## Optimization gadget
    palindromic_windows = []

    covered_indices = set()

    for center in range(len(car_path)):
        left_bound = center
        right_bound = center
        dist = 1
        while (center-dist >= 0) and (center+dist < len(car_path)) and (car_path[center-dist] == car_path[center+dist]):
                left_bound -= 1
                right_bound += 1
                dist += 1
        contained = True
        for i in range(left_bound, right_bound+1):
            if i not in covered_indices:
                contained = False
        if not contained:
            palindromic_windows.append((left_bound, right_bound))
            for i in range(left_bound, right_bound+1):
                covered_indices.add(i)

    clean_windows = []

    for window1 in palindromic_windows:
        contained = False
        for window2 in palindromic_windows:
            if window1 != window2:
                if (window1[0] >= window2[0]) and (window1[1] <= window2[1]):
                    contained = True
                    break
        if not contained:
            clean_windows.append(window1)

    clean_path = [car_path[clean_windows[0][0]]]

    home_set = set(home_list)

    for window in clean_windows:
        if clean_path[-1] != car_path[window[0]]:
            clean_path.append(car_path[window[0]])
        for i in range(window[0] + 1 + ((window[1]-window[0])//2), window[1]):
            if car_path[i] in home_set:
                #print(i)
                dropoff_mappings[car_path[i]] = [car_path[window[0] + ((window[1]-window[0])//2)]]
                for j in range(window[1]-1, window[0] + 1 + ((window[1]-window[0])//2), -1):
                    clean_path.append(car_path[j])
                for j in range(window[0] + 1 + ((window[1]-window[0])//2), window[1]):
                    clean_path.append(car_path[j])
                clean_path.append(car_path[window[0]])
                break

    for home in clean_path:
        if home in home_list:
            if home in dropoff_mappings:
                dropoff_mappings[home].append(home)
                home_list.remove(home)
            else: 
                dropoff_mappings[home] = [home]
                home_list.remove(home)
                
    node_list = clean_path

    return node_list, dropoff_mappings

def trash_solve(list_of_locations, list_of_homes, starting_car_location, adjacency_matrix, params=[]):
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
    G = nx.DiGraph()

    start_node = list_of_locations.index(starting_car_location)

    for i in range(num_locations):
        for j in range(num_locations):
            if (adjacency_matrix[i][j] != 'x'):
                G.add_edge(i, j, weight=int(adjacency_matrix[i][j]))

    mod_G = make_gadget_graph(G)

    E = G.edges
    mod_E = mod_G.edges

    paths, sp_distances = nx.floyd_warshall_predecessor_and_distance(G)

    

    # Initializing capacities
    edge_capacity_dict = {}
    for edge in mod_E:
        edge_capacity_dict[edge] = 1 
    edge_capacity_dict[('s', start_node)] = float('inf')

    m = Model()


    ################## VARIABLES ##################

    # Edge inclusion variables.
    # Edge (i,j) in subtour => edge_vars[i,j] = 1
    edge_vars = {}

    for edge in E:
        edge_vars[edge] = m.addVar(name="e_"+str(edge), vtype=GRB.BINARY)

    m.update()
    
    # TA dropoff variables.
    # TA t dropped off at location l => ta_vars[t,l] = 1
    ta_vars = {}

    for t in range(num_homes):
        for l in range(num_locations):
            ta_vars[t,l] = m.addVar(name="ta_"+str(t)+"@"+str(list_of_locations[l]), vtype=GRB.BINARY)
    
    m.update()
    
    # Edge flow variables.
    # Edge (i, j) in tour that includes start => flow_vars[i,j] = 1
    flow_vars = {}

    for edge in mod_E:
        flow_vars[edge] = m.addVar(name="flow_"+str(edge), vtype=GRB.BINARY)
    flow_vars[('s', start_node)] = m.addVar(name="flow_source", vtype=GRB.INTEGER)
    
    for node in mod_G.nodes:
        if type(node) != int and "mid" in node:
            flow_vars[(node, 't')] = m.addVar(name="flow_target", vtype=GRB.INTEGER)

    m.update()
    
    ################# CONSTRAINTS #################

    # Cycle enforcement constraint
    for v in range(num_locations):
        m.addConstr(quicksum(edge_vars[(i,v)] if (i,v) in E else 0 for i in G.nodes) == quicksum(edge_vars[(v,j)] if (v,j) in E else 0 for j in G.nodes), "indegree_equals_outdegree")
    
    m.update()

    # One dropoff per TA constraint
    for t in range(num_homes):
        m.addConstr(quicksum(ta_vars[t,l] for l in range(num_locations)) == 1, "one_dropoff_ta" + str(t))
    
    m.update()

    # All TAs dropped off
    m.addConstr(quicksum(ta_vars.values()) == num_homes, "all_dropped_off")

    m.update()

    # TA dropped off on path
    for t in range(num_homes):
        for l in range(num_locations):
            m.addConstr(ta_vars[t,l] <= quicksum([edge_vars[l, v] if (l, v) in E else 0 for v in range(num_locations)]), "ta" + str(t) + "_loc" + str(l) + "_valid_dropoff")

    m.update()

    # Source flow constraint
    m.addConstr(flow_vars[('s', start_node)] == 1, "source_flow")

    # Flow i,j is less than Capacity i,j
    for edge in mod_E:
        m.addConstr(flow_vars[edge] <= edge_capacity_dict[edge], "capacity_not_exceeded")

    m.update()

    # Flow passed along path from node to node
    for v in mod_G.nodes:
        if v != start_node:
            m.addConstr(quicksum(flow_vars[(i,v)] if (i,v) in mod_E else 0 for i in mod_G.nodes) == quicksum(flow_vars[(v,j)] if (v,j) in mod_E else 0 for j in mod_G.nodes), "flow_conserved")
        if type(v) == int:
            m.addConstr(flow_vars[(str(v)+"_mid_L", 't')] + flow_vars[(str(v)+"_mid_R", 't')] == 1)

    m.addConstr(flow_vars[('s', start_node)] + quicksum(flow_vars[(i,start_node)] if (i,start_node) in mod_E else 0 for i in mod_G.nodes)  == quicksum(flow_vars[(start_node, v)] if (start_node, v) in mod_E else 0 for v in mod_G.nodes), "flow_start")

    m.update()

    # Flow is non-negative
    for flow in flow_vars.values():
        m.addConstr(flow >= 0, "non-negative flow")
    
    m.update()

    # Edge included => Flow == 1 == Inclusion Variable
    for edge in E:
        m.addConstr(flow_vars[(str(edge[0])+"_out", edge[1])] == edge_vars[edge], "inclusion_implies_flow")

    m.update()

    ################Setting Objective###############
    
    m.setObjective((2.0/3.0)*quicksum(edge_vars[e]*int(adjacency_matrix[e[0]][e[1]]) for e in E) + quicksum(quicksum(ta_vars[t,l]*sp_distances[t][l] for l in range(num_locations)) for t in range(num_homes)), GRB.MINIMIZE)
    
    m.update()

    m.optimize()

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
