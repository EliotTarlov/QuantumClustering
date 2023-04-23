import localsolver
import sys
from QCutting import distance, printClustering, createGraph
from networkx import get_edge_attributes
def generate_graph(points):
    G=createGraph(points)
    print(get_edge_attributes(G,'weight'))
    return len(G.edges),len(G.nodes),[x[0] for x in G.edges],[x[1] for x in G.edges], get_edge_attributes(G,'weight').values()
def main(instance_file, output_file, time_limit):
    from ast import literal_eval
    with open(instance_file) as f:
         points=literal_eval(f.read())
    n, m, origin, dest, w = generate_graph(points)
    with localsolver.LocalSolver() as ls:
        # Declare the optimization model
        model = ls.model
        # Decision variables x[i]
        # True if vertex x[i] is on the right side of the cut
        # and false if it is on the left side of the cut
        x = [model.bool() for i in range(n)]
        # An edge is in the cut-set if it has an extremity in each class of the bipartition
        incut = [None] * m
        for e in range(m):
            incut[e] = model.neq(x[origin[e] - 1], x[dest[e] - 1])
        # Size of the cut
        cut_weight = model.sum(w[e] * incut[e] for e in range(m))
        model.maximize(cut_weight)
        model.close()
        # Parameterize the solver
        ls.param.time_limit = time_limit
        ls.solve()
        # Write the solution in a file with the following format:
        #  - objective value
        #  - each line contains a vertex number and its subset (1 for S, 0 for V-S)
        if output_file != None:
            with open(output_file, 'w') as f:
                f.write("%d\n" % cut_weight.value)
                # Note: in the instances the indices start at 1
                for i in range(n):
                    f.write("%d %d\n" % (i + 1, x[i].value))
        else:
            from QCutting.py import printClustering

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python maxcut.py inputFile [outputFile] [timeLimit]")
        sys.exit(1)

    instance_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) >= 3 else None
    time_limit = int(sys.argv[3]) if len(sys.argv) >= 4 else 10
    main(instance_file, output_file, time_limit)
