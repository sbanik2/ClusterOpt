with open("dumpfile.dat", "r") as infile, open("plot.mcts", "w") as outfile:
#with open("mae.dat", "r") as infile, open("plot.mcts", "w") as outfile:
    minval = 1e7
    coords = ""
    cnt = 0
    minloc = 0
    for i, line in enumerate(infile):
        cnt += 1
        col = line.split("|")
#        col = line.split("|")
        eng = float(col[-1])
        if minval > eng:
            minval = eng
            minloc = i+1
            coords = col[1]
            outfile.write("%s %s\n"%(i+1, minval))

print("Lowest Value: %s"%(minval))
#print("Lowest Coords: %s"%(coords))
print("Number of Total Evaluations: %s"%(cnt))
print("Number of Evaluations till Minima was found: %s"%(minloc))
