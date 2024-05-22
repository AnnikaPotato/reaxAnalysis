import re
import os.path
from collections import deque as deque
from random import randrange

class FileParameters():
    def __init__(self, sourcefile = '../parameters.txt') -> None:
        para = dict()
        with open(sourcefile) as fhand:
            for line in fhand:
                if not line.strip(): continue
                line = line.split('=')
                para_key = line[0].strip().lower()
                if para_key not in para:
                    para[para_key] = [line[-1].strip()]
                else:
                    para[para_key].append(line[-1].strip())

        with open(para['in_file'][0],'r') as fhand:
            for line in fhand:
                if not line.lstrip().startswith('read_data'):
                    continue
                line = re.sub('\s{2,}', ' ', line)
                line = re.sub('\t', ' ', line)
                sourcefilename = line.rstrip().split(' ')[-1]
                para['data_file'] = sourcefilename
                pos = sourcefilename.find('.')
                para['curr_run'] = re.match('[a-zA-Z]+', sourcefilename[:pos+1])[0]
                curr_cycle = re.search('[0-9]+', sourcefilename[:pos+1])
                if curr_cycle:
                    para['curr_cycle'] = int(curr_cycle[0])
                else:
                    para['curr_cycle'] = 1
                break

        self.para = para
        #return self.para

class Atom():
    """
    Analyze the input data file, catagorize the atom and save the indices of the same atom in a file
    """
    def __init__(self, fname: str) -> None:
        self.fname = fname #save the name of the data file
        self.atoms = []
        self.names = []
        self.corresponding = []

    def extract(self, targetfolder = './atom/') -> int:
        if not os.path.isdir(targetfolder):
            os.mkdir(targetfolder) #create target folder if not exists

        with open(self.fname, "r") as fhand:
            ## read the file line by line and find the number of atom types
            for line in fhand:
                num_atom_type = re.findall("([0-9]+) atom types", line)
                if num_atom_type: break
            if not num_atom_type:
                print("Unable to find the number of atom types. Check the input file")
                return 1
            num_atom_type = int(num_atom_type[0])

            ## read the file line by line and find the name of each atom type
            for line in fhand:
                if line.strip().startswith("Masses"): break
            findname = False

            i = 0
            for line in fhand:
                atomname = re.findall("#\s*([0-9a-zA-Z]+)", line)
                if atomname:
                    self.names.append(atomname[0])
                    i += 1
                if i == num_atom_type:
                    findname = True
                    break
            for i in range(len(self.names)):
                self.atoms.append([])
            ## continue to read the file and catagorize each atom
            findstart = False
            for line in fhand:
                line = line.strip()
                if not findstart and not line.startswith('1'): continue
                findstart = True
                line = line.split()
                index = int(line[1]) - 1
                if  index < num_atom_type:
                    self.atoms[index].append(str(int(line[0]) - 1))
                    self.corresponding.append(index)
                else:
                    print("The number of atoms claims at the beginning of the file doesn't match.")
                    return 1

        ## print the catagorized atom indices into separated files
        for i in range(num_atom_type):
            if findname:
                new_fname = targetfolder + self.names[i]
            else:
                new_fname = targetfolder + '#' + str(i+1)
            contents = " ".join(self.atoms[i])
            with open(new_fname, 'w') as fhand:
                fhand.write(contents)
        return 0

class ExtractTrajectory():
    """
    separate and save trajectory files timestepwisely
    """
    def __init__(self, file = 'trajectory.lammpstrj', targetfolder = './trj/') -> None:
        if not os.path.isdir(targetfolder):
            os.mkdir(targetfolder)

        contents = []
        FindHeader = False
        with open(file, 'r') as fhand:
            for line in fhand:
                line = line.strip()
                if line.startswith('ITEM: TIMESTEP'):
                    FindHeader = True
                    if contents:
                        contents.sort()
                        for _, v in contents:
                            new_file_hand.write(v)
                        new_file_hand.close()
                        contents = []
                    continue

                if FindHeader:
                    FindHeader = False
                    timestep = line.strip().split(' ')
                    new_file = targetfolder + timestep[0] +'.lammpstrj'
                    new_file_hand = open(new_file, 'w')
                    new_file_hand.write('ITEM: TIMESTEP \n' + timestep[0] + '\n')
                    continue

                res = re.findall('^[0-9]+ [\-0-9]+.',line)
                if not res:
                    new_file_hand.write(line +'\n')
                else:
                    key = int(line.split(' ')[0])
                    contents.append((key, line + '\n'))
        contents.sort()
        for _, v in contents:
            new_file_hand.write(v)
        new_file_hand.close()

class ExtractBond():
    """
    separate and save bond information files timestepwisely
    """
    def __init__(self, file = 'bonds.reaxc', targetfolder = './bond/') -> None:
        if not os.path.isdir(targetfolder):
            os.mkdir(targetfolder)
        contents = []

        with open(file,'r') as fhand:
            for line in fhand:
                line = line.strip()

                if line.startswith('# Timestep'):
                    line = line.split(' ')
                    if contents:
                        contents.sort()
                        new_file = targetfolder + 'bond.ts.' + timestep
                        with open(new_file, 'w') as new_file_hand:
                            for _, content in contents:
                                new_file_hand.write(content +'\n')
                        contents = []
                    timestep = line[2]
                    continue

                if line.startswith('#'): continue

                line = re.sub('\s{2,}', ' ', line)                  #remove duplicated spaces in the line
                atom_index = int(line.split(' ')[0])
                contents.append((atom_index, line))

            new_file = targetfolder + 'bond.ts.' + timestep
            contents.sort()
            with open(new_file, 'w') as new_file_hand:
                for _, content in contents:
                    new_file_hand.write(content +'\n')
        return

class AnalyzeBond():
    """
    Analyze the extracted bond files made by ExtractBond
    Get and save all connections between two atoms
    The saved connection information can be used in VMD with topotools
    Note that in VMD, atom index starts from 0 while in LAMMPS from 1
    """
    def __init__(self, source_folder = './bond/', target_folder = './connections/', trj_folder = './trj/') -> None:
        self.source_folder = source_folder
        self.target_folder = target_folder
        self.trj_folder = trj_folder

    def connections(self, startstep = 0, endstep = 0, stepwidth = 1, cutoff = 0) -> None:
        if not os.path.isdir(self.target_folder):
            os.mkdir(self.target_folder)
        if cutoff and not os.path.isdir('./conn_ncb/'):
            os.mkdir('./conn_ncb/')
        if cutoff:
            cutoff = cutoff ** 2

        for step in range(startstep, endstep + 1, stepwidth):
            step = str(step)
            fname = self.source_folder + 'bond.ts.' + step
            num_atom = 0
            contents = []
            with open(fname, 'r') as fhand:
                for line in fhand:
                    num_atom += 1
                    line = line.strip().split(' ')
                    num_nei = int(line[2]) # number of other atoms connected to the atom
                    atom_index = int(line[0]) - 1
                    for i in range(num_nei):
                        nei_index = int(line[i + 3]) - 1
                        if nei_index > atom_index:
                            contents.append((nei_index, atom_index))
            new_file = self.target_folder + step
            with open(new_file, 'w') as new_file_hand:
                for atom1, atom2 in contents:
                    new_file_hand.write('topo addbond ' + str(atom1) + ' ' + str(atom2) + '\n')

            #if trajectory files extracted from ExtractTrajectory are not found, bond information with cutoff
            #will be dismissed
            #'bond information with cutoff' means that if the bond length bewteen two atoms is longer than
            #the given cutoff value, it will not be print
            #suitable for visualization in VMD
            #cutoff = 0 means unable the function

            if not cutoff: continue
            fname = self.trj_folder + step + '.lammpstrj'
            try:
                fhand = open(fname)
            except:
                cutoff = 0
                print('''trajectory files extracted from ExtractTrajectory are not found,
                bond information with cutoff is dismissed''')
                continue
            fhand.close()

            atom_coordinates = [(0,0,0)] * num_atom
            findstart = False
            with open(fname, 'r') as fhand:
                for line in fhand:
                    if line.startswith('ITEM: ATOMS'):
                        findstart = True
                        continue
                    if not findstart: continue
                    line = line.split(' ')
                    if len(line) < 4: continue
                    atom_index = int(line[0]) - 1
                    atom_coordinates[atom_index] = (float(line[1]), float(line[2]), float(line[3]))

            new_file = './conn_ncb/' + step
            with open(new_file, 'w') as new_file_hand:
                for atom1, atom2 in contents:
                    dist = 0
                    for i in range(3):
                        dist += (atom_coordinates[atom1][i] - atom_coordinates[atom2][i]) ** 2
                    if dist >= cutoff: continue
                    new_file_hand.write('topo addbond ' + str(atom1) + ' ' + str(atom2) + '\n')

class Fragments():
    def __init__(self, targetstep: str, write = 1) -> None:
        fname = 'bonds.reaxc'
        try:
            fhand = open(fname,'r')
        except:
            print('Cannot open the source file in ' + fname +'. Please check.')
            return
        i = 0
        connections = dict()
        find = False
        header = '# Timestep ' + targetstep
        for line in fhand:
            if not find and line.startswith('#'):
                if line.startswith(header):
                    find = True
                    continue
                else:
                    continue
            if not find: continue
            if line.lstrip().startswith('#'): continue
            line = line.lstrip().split(' ')
            counts = int(line[2])
            for index in range(counts):
                atom = str(int(line[0]) - 1)
                nei = str(int(line[index + 3]) - 1)
                if atom not in connections:
                    connections[atom] = [nei]
                else:
                    connections[atom].append(nei)
        fhand.close()
        self.connections = connections

        stack = deque()
        seen = set()
        fragments = []
        for atom in range(len(connections)):
            atom = str(atom)
            if atom in seen: continue
            fragment = []
            fragment.append(atom)
            seen.add(atom)
            stack.append(atom)
            while stack:
                atom = stack.popleft()
                if atom not in connections: continue
                connections[atom].sort()
                for nei in connections[atom]:
                    if nei in seen: continue
                    stack.append(nei)
                    fragment.append(nei)
                    seen.add(nei)
            fragments.append(fragment[:])
        self.fragments = fragments

        if not write: return
        with open('all_fragments.txt','w') as new_file_hand:
            for i in range(len(fragments)):
                new_file_hand.write('Fragment ' + str(i) +':\n')
                new_file_hand.write(' '.join(fragments[i]))
                new_file_hand.write('\n')


    def DeleteSmallMolecules(self, cutoff: int, curr_run: str, next_run: str, end_step: str, namelist: list, names: str) -> None:
        def print_delete_fragments(j):
            res = []
            counts = [0] * len(names)
            for atom in self.fragments[j]:
                counts[namelist[int(atom)]] += 1
            for i in range(len(counts)):
                if not counts[i]: continue
                res.append(names[i])
                if counts[i] > 1:
                    res.append(str(counts[i]))
            return ''.join(res)

        to_delete = set()
        with open('deleted_fragments.txt','w') as fhand:
            for i in range(len(self.fragments)):
                if len(self.fragments[i]) > cutoff:
                    #print(i)
                    continue
                for atom in self.fragments[i]:
                    to_delete.add(atom)
                to_write = print_delete_fragments(i) + '\n'
                fhand.write(to_write)
        source = './trj/' + end_step +'.lammpstrj'
        output = '../' + next_run +'/' + next_run + '.data'
        path = '../' +next_run + '/'
        if not os.path.isdir(path):
            os.mkdir(path)
        source_hand = open(source, 'r')
        data_file= curr_run + '.data'
        find = False
        count = 1
        suffix = ["xlo xhi\n", "ylo yhi\n", "zlo zhi\n"]

        with open(output, 'w') as new_file_hand:
            new_file_hand.write("#\n")
            for line in source_hand:
                if not find:
                    num_atm = re.findall("ITEM: NUMBER OF ATOMS", line)
                    if num_atm:
                        break
            num_atm = int(source_hand.readline()) - len(to_delete)
            new_file_hand.write( str(num_atm)+ " atoms\n")
            num_atm = str(len(names))
            new_file_hand.write(num_atm + " atom types\n")
            source_hand.readline()

            for i in range(3):
                temp = source_hand.readline().strip()
                temp = temp + " " + suffix[i]
                new_file_hand.write(temp)

            new_file_hand.write("\nMasses\n\n")

            with open(data_file, 'r') as data_hand:
                for line in data_hand:
                    if not line.startswith("Masses"):
                        continue
                    else:
                        break
                data_hand.readline()
                for line in data_hand:
                    if line != "\n":
                        new_file_hand.write(line)
                    else:
                        break

            new_file_hand.write("\nAtoms\t# charge\n\n")

            source_hand.readline()
            for line in source_hand:
                line = line.lstrip().split(' ')
                if str(int(line[0]) - 1) in to_delete: continue
                to_write = str(count)
                count += 1
                to_write += " " + str(namelist[int(line[0]) - 1] + 1) + " 0 "
                to_write += " ".join(line[1:])
                new_file_hand.write(to_write)
        source_hand.close()


class AtomCount():
    def __init__(self, cutoff: int, atomfile: str, fragmentfile = 'all_fragments.txt', ) -> None:
        total = set()
        with open(fragmentfile,'r') as fhand:
            for line in fhand:
                if line.startswith('Fragment'): continue
                line = line.split()
                if len(line) <= cutoff: continue
                for item in line:
                    total.add(str(int(item)+1))
        with open(atomfile, 'r') as fhand:
            for line in fhand:
                num_atom_type = re.findall("([0-9]+) atom types", line)
                if num_atom_type: break
            num_atom_type = int(num_atom_type[0])
            counts = dict()
            names = []
            i = 0
            for line in fhand:
                atomname = re.findall("#\s*([0-9a-zA-Z]+)", line)
                if atomname:
                    counts[str(i+1)] = 0
                    names.append(atomname[0])
                    i += 1
                if i == num_atom_type:
                    break
            findstart = False
            for line in fhand:
                line = line.strip()
                if not findstart and not line.startswith('1'): continue
                findstart = True
                line = line.split()
                if line[0] in total:
                    counts[line[1]] += 1
        for i in range(num_atom_type):
            print(names[i], end = '')
            print(counts[str(i+1)], end = '')
        print('')


class NewInfile():
    def __init__(self) -> None:
        pass
    def create(self, next_run:str, temp_incr: str, end_step:str) -> None:
        temp_incr = int(temp_incr)
        curr_hand = open('in.lammps','r')
        path = '../' + next_run +'/'+ 'in.lammps'
        next_hand = open(path, 'w')

        for line in curr_hand:
            if line.lstrip().startswith('read_data'):
                    next_hand.write('read_data\t'+ next_run +'.data\n')
                    break
            else:
                next_hand.write(line)


        for line in curr_hand:
            if line.lstrip().startswith('velocity'):
                last_temp = self.temp_read(end_step)
                seed = str(randrange(10000000))
                next_hand.write('velocity all create ' + last_temp + ' ' + seed + ' rot yes dist gaussian\n')
                break
            else:
                next_hand.write(line)

        for line in curr_hand:
            if re.findall('npt',line):
                line = re.sub('\s{2,}', ' ', line)
                line = line.split(' ')
                line[5] = str(int(line[5]) + temp_incr)
                line[6] = str(int(line[6]) + temp_incr)

                next_hand.write(' '.join(line))
                break
            else:
                next_hand.write(line)

        for line in curr_hand:
            next_hand.write(line)



        next_hand.close()
        curr_hand.close()

    def temp_read(self, end_step: str) -> str:

        with open('log.lammps','r') as fhand:
            for line in fhand:
                if re.findall('[0-9]+\s+[.0-9\-]+\s+[.0-9\-]+\s+[.0-9\-]+\s+[.0-9\-]+', line):
                    line = re.sub('\s{2,}', ' ', line)
                    line = line.lstrip().split(' ')
                    if line[0] == end_step:
                        return line[1]
