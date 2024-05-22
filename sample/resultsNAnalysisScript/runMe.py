import reax_ana

info = reax_ana.FileParameters('./parameters.txt')
start_step = int(info.para['start_step'][0])
end_step = int(info.para['end_step'][0])
stepwidth = int(info.para['stepwidth'][0])
cutoff = int(info.para['cutoff'][0])

for i in range(len(info.para['runnames'])):
    if info.para['runnames'][i] == info.para['curr_run']:
        curr_run_index = i
        break
if info.para['curr_cycle'] == int(info.para['cycles'][i]):
    if curr_run_index < (len(info.para['runnames']) - 1):
        next_run_index = curr_run_index + 1
        next_cycle_num = 1
    else: curr_run_index = -1
else:
    next_run_index = curr_run_index
    next_cycle_num = info.para['curr_cycle'] + 1

if curr_run_index >= 0:
    next_cycle = info.para['runnames'][next_run_index] + str(next_cycle_num)
    curr_cycle = info.para['runnames'][curr_run_index] + str(info.para['curr_cycle'])


atoms = reax_ana.Atom(info.para['data_file'])
atoms.extract()
reax_ana.ExtractTrajectory()
reax_ana.ExtractBond()
reax_ana.AnalyzeBond().connections(start_step, end_step, stepwidth, cutoff)

if curr_run_index >= 0:
    reax_ana.Fragments(str(end_step), 1).DeleteSmallMolecules(cutoff, curr_cycle, next_cycle, str(end_step), atoms.corresponding, atoms.names)
    reax_ana.NewInfile().create(next_cycle, info.para['temp_incr'][next_run_index], str(end_step))
else:
    reax_ana.Fragments(str(end_step), 1)
