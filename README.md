Efficient Public Transport Planning on Roads
========================================================================

This repository contains the source code and data used in the paper. 

Usage
---------------

### Environment

gcc version: 9.4.0 

OS: Ubuntu/CentOS (Unix)

### Compilation

make

### Execution

We directly use the results of vk-TSP and ETA-Pre stored in the folder KResults and QResults for varying K and Q, respectively. Our algorithm EBRR could be executed by the following statement.

./EBR [network path] [transit path] [query path] [$\alpha$] [K] [D] [path of ETA-Pre results (optional)] [path of vk-TSP results (optional)]

The default setting uses

network path=data/Chicago/network.txt

transit path=data/Chicago/transit.txt

query path=data/Chicago/queries.txt

$\alpha$, K, D=2000, 30, 2

path of ETA-Pre results (optional)=KResults/ChiResults/eta-pre/30/route_stops.txt

path of vk-TSP results (optional)=KResults/ChiResults/vk-tsp/30/route_stops.txt

A complete statement is ./EBR data/Chicago/network.txt data/Chicago/transit.txt data/Chicago/queries.txt 2000 30 2 KResults/ChiResults/eta-pre/30/route_stops.txt KResults/ChiResults/vk-tsp/30/route_stops.txt

The final results are stored in a file named "results.txt" under the same folder as the code. The results.txt gives the number of stops $|B_{r*}|$ in the new bus route in the first line. Then, the results of EBRR, ETA-Pre, and Vk-TSP are given in the next starting with EBRR, ETA-Pre, and Vk-TSP, respectively. Finally, there is a new line showing the running time of EBRR.

When we run the program, the screen will show some additional information in different steps of the EBRR, such as the preprocessing, the selection, and the final pruning. 

### Data Description

The data used in the experiments are stored in the data folder. There are three folders under the data folder corresponding to the three cities: Chicago, NYC, and Orlando. Each folder contains its network, transit, and queries files needed by the algorithms. Besides, since we consider different query sets in different local areas of the city, we include the corresponding data in the sub-folders named 1-4. 

Each network file contains the number of nodes $n$ and edges $m$ in the first line. The following $n$ lines give the longitude and latitude of the nodes. The following m lines have the origin, the destination node, and the weights of edges w_e. Note that one unit of edge weights is one kilometer.

Each transit file contains the number of new and existing stops in the first line. The following lines first show information of each new stop and then each existing stop. For each new stop in each line, we give the edge id, the distance to the origin node of the edge, the longitude and latitude of each new stop. For each existing stop in each line, the same information as a new stop is first given and followed by the set $routes(v)$ of existing routes that pass through the existing stop. 

Each query file contains the number of queries |Q| in the first line. The following |Q| lines gives the corresponding query nodes.
