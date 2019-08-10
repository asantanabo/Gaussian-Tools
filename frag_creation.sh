import sys

file_path = sys.argv[1] #Strips filename from terminal
file= open(file_path)
contents = file.readlines() #Gets coordinate data from .xyz file
file.close()
i = 0
atom = []
for lines in contents: #Arranges data into name/x/y/z/number
    if i > 1:
        atom = atom + [[contents[i].split()[0], float(contents[i].split()[1]), float(contents[i].split()[2]), float(contents[i].split()[3]), i-2]]
    i = i + 1

vdw_radii = {
    "H":  0.33,
    "C":  0.76,
    "O":  0.66,
    "N":  0.71,
    "Pd": 1.39,
    "PD": 1.39,
    "F":  0.61, #mod
    "S": 1.10   #mod
    }

fragments = []
def bonding(point, fragments):
   #Scans all other atoms not in current fragment list to see if within bonding distance of selected atomic_mass
   #Adds all newly found bonding atoms to the current fragment list
    test = point
    for x in atom:
        if test[4] != x[4]:
            sep_dist = ((test[1]-x[1])**2 + (test[2]-x[2])**2 + (test[3]-x[3])**2)**0.5
            if sep_dist <= (vdw_radii[test[0]] + vdw_radii[x[0]]):
                if x[4] not in fragments[-1]:
                    fragments[-1] += [x[4]]
    return point

total_length = 0
Found = False
for x in fragments:
    total_length += len(x)

while total_length != len(atom): #ie until all atoms are asssigned to a list

    for x in atom:
        Found = False
        for y in fragments: #Checks if current atom has been assigned to a list
            if x[4] in y:
                Found = True
        if Found == False: #Works on first atom not found in a list
            next = x[4]
            fragments += [[next]] #Adds this atom to a new list, as doesn't belong to any other
            break
    print(next)
    for x in fragments[-1]: #Works on the most recently created fragement list
        bonding(atom[x], fragments) #Will iterate until fragement list stops being appended to, i.e. when completed
    total_length = 0
    for x in fragments:
        total_length += len(x) #Counts assigned atoms

print(fragments)

i = 0
for y in fragments: #Prints each fragment to its own .xyz file
    i += 1
    f = open(file_path[:-4] + '_frag' + str(i) + '.xyz',"w+")
    f.write(str(len(y)) + '\n')
    f.write('\n')
    for x in y:
        f.write(atom[x][0] + ' ' + str(atom[x][1]) + ' ' + str(atom[x][2]) + ' ' + str(atom[x][3]) + '\n')
    f.close()
