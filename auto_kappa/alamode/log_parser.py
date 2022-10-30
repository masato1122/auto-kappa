# -*- coding: utf-8 -*-
import os.path
import numpy as np
from optparse import OptionParser
import ase.io
import datetime
import yaml
import warnings
import glob

from auto_kappa import output_directories as out_dirs

class AkLog():
    """
    
    How to Use
    -------------
    >>> directory = './mp-149'
    >>> log = AkLog(directory)
    >>>
    >>> ## file names are automatically set.
    >>> ## "log.yaml" and "fig_times.png" are created in directory+'/result'.
    >>> log.write_yaml() 
    >>> log.plot_times()
    >>>
    >>> log.write_yaml(outfile='log.yaml')
    >>> log.plot_times(figname='fig_times.png')
    
    """
    def __init__(self, directory):
        
        self._directory = directory
        self._out = get_ak_logs(self._directory)
    
    @property
    def directory(self):
        return self._directory
    
    @property
    def out(self):
        return self._out
    
    def write_yaml(self, outfile=None):
        if outfile is None:
            outfile = self.directory+'/'+out_dirs['result']+'/log.yaml'
        with open(outfile, "w") as f:
            yaml.dump(self.out, f)
            print("")
            print(" Output", outfile)
    
    def get_times(self):
        
        times = {}
        for k1 in self.out:

            if isinstance(self.out[k1], dict) == False:
                continue
            
            if 'time' in self.out[k1]:
                times[k1] = self.out[k1]['time']['value']
            
            for k2 in self.out[k1]:

                if isinstance(self.out[k1][k2], dict) == False:
                    continue
            
                if 'time' in self.out[k1][k2]:
                    lab = k1+'_'+k2
                    times[lab] = self.out[k1][k2]['time']['value']
        
        return times

    def plot_times(self, figname=None):

        if figname is None:
            figname = self.directory+'/'+out_dirs['result']+'/fig_times.png'
        
        nonskip_labels = ['force', 'cv(lasso)', 'lasso(lasso)', 'kappa']
        
        all_times = self.get_times()

        times = []
        labels = []
        time_others = 0.
        for key in all_times:
             
            if '_' in key:
                data = key.split('_')
                if data[0] == 'kappa':
                    lab = "%s(%s)" % (data[0], data[1])
                else:
                    lab = "%s(%s)" % (data[1], data[0])
            else:
                lab = key
            
            flag_skip = True
            for ll in nonskip_labels:
                if ll in lab:
                    flag_skip = False
            
            if flag_skip:
                time_others += all_times[key]
            else:
                times.append(all_times[key])
                labels.append(lab)
        
        times.append(time_others)
        labels.append('others')

        from auto_kappa.plot.pltalm import plot_times_with_pie
        plot_times_with_pie(times, labels, figname=figname)
        

def _extract_data(filename, word, back_id=-1):
    """ Return data for ``word`` in ``filename``
    """
    lines = open(filename, 'r').readlines()
    values_get = []
    for line in lines:
        if word.lower() in line.lower():
            data = line.split()[back_id]
            values_get.append(float(data))
    if len(values_get) == 0:
        return None
    else:
        return values_get

def _extract_lines(filename, word):
    """ Return lines for ``word`` in ``filename``
    """
    lines = open(filename, 'r').readlines()
    lines_get = []
    for line in lines:
        if word.lower() in line.lower():
            lines_get.append(line.replace("\n", ""))
    if len(lines_get) == 0:
        return None
    else:
        return lines_get

#class Value:
#    def __init__(self, value, unit):
#        self._value = value
#        self._unit = unit
#    @property
#    def value(self):
#        return self._value
#    @property
#    def unit(self):
#        return self._unit
#    def __repr__(self):
#        line = "%f %s" % (self.value, self.unit)
#        return line

def _get_alamode_runtime(filename):
    
    def _adjust_time_line(line_orig):
        data = line_orig.split()
        line_mod = ""
        for i, dd in enumerate(data):
            if i == 2:
                line_mod += "%02d" % int(dd)
            else:
                line_mod += dd
            if i != len(data)-1:
                line_mod += " "
        return line_mod
    try:
        line0 = _extract_lines(
                filename, 'job started at')[0].split('started at')[1].replace('\n', "")
        line1 = _extract_lines(
                filename, 'job finished at')[-1].split('finished at')[1].replace("\n", "")
        
        time0 = datetime.datetime.strptime(
                _adjust_time_line(line0), "%a %b %d %H:%M:%S %Y")
        time1 = datetime.datetime.strptime(
                _adjust_time_line(line1), "%a %b %d %H:%M:%S %Y")
        diff = time1 - time0
        runtime = {'value': diff.seconds, 'unit': 'sec'}
        return runtime
    
    except Exception:
        return None

def read_log_fc(filename):
    
    if os.path.exists(filename) == False:
        return None
    
    out = {}
    v = _get_alamode_runtime(filename)
    if v is not None:
        out['time'] = v
     
    lines = open(filename, 'r').readlines()
    
    ### read different parameters
    i0_cut = None
    for ii, line in enumerate(lines):

        if "NKD =" in line:
            data = line.split()
            for j, dd in enumerate(data):
                if dd == "NKD":
                    out['nkd'] = int(data[j+2])
        
        if "NORDER =" in line:
            data = line.split()
            out['norder'] = int(data[-1])
        
        if "NBODY =" in line:
            data = line.split()
            out['nbody'] = []
            for i in range(out['norder']):
                out['nbody'].append(int(data[i+2]))
        
        if "cutoff radii matrix" in line.lower():
            i0_cut = ii
    
        if "Atomic species:" in line:
            i0_species = ii
    
    ## get atomic species
    try:
        if i0_species is not None:
            species = []
            for ik in range(out['nkd']):
                data = lines[i0_species+1+ik].split()
                species.append(data[1])
        out['atomic_species'] = species
    except Exception:
        out['atomic_species'] = None
    
    ## get cutoff radii
    try:
        cutoff_mat = {'unit': 'Bohr', 'value': []}
        for iorder in range(out['norder']):
            cutoff_mat['value'].append([])
            for i1 in range(out['nkd']):
                num = i0_cut + iorder * (out['nkd'] + 2) + i1 + 2
                data = lines[num].split()
                cutoff_mat['value'][iorder].append([])
                for dd in data:
                    if dd.lower() == 'none':
                        cutoff_mat['value'][iorder][i1].append(None)
                    else:
                        cutoff_mat['value'][iorder][i1].append(float(dd))
        out['cutoff_radii_matrix'] = cutoff_mat
    except Exception:
        pass
    
    ## fitting error
    v = _extract_data(filename, 'fitting error (%)')
    if v is not None:
        out['fitting_error'] = {'value': float(v[-1]), 'unit': "%"}
    
    ## residual
    v = _extract_data(filename, "residual (%)")
    if v is not None:
        out['residual'] = {'value': float(v[0]), 'unit': "%"}
    
    ## number of structures
    v = _extract_data(
            filename, "entries will be used for training", back_id=-7)
    if v is not None:
        out['number_of_structures'] = int(v[0])
     
    ## number of parameters
    v = _extract_data(filename, "Total Number of Parameters")
    if v is not None:
        out['number_of_parameters'] = int(v[0])
     
    ## number of free parameters
    v = _extract_data(filename, "Total Number of Free Parameters")
    if v is not None:
        out['number_of_free_parameters'] = int(v[0])
     
    ### number of FCs
    vs = _extract_data(filename, "Number of  HARMONIC FCs")
    if vs is not None:
        out['number_of_fc2'] = int(vs[0])
    #
    for order in range(3, 10):
        vs = _extract_data(filename, "Number of   ANHARM%d FCs" % order)
        if vs is not None:
            out['number_of_fc%d' % order] = int(vs[0])
    
    ### Number of free FCs
    vs = _extract_data(filename, "Number of free HARMONIC FCs")
    if vs is not None:
        out['number_of_free_fc2'] = int(vs[0])
    #
    for order in range(3, 10):
        vs = _extract_data(filename, "Number of free  ANHARM%d FCs" % order)
        if vs is not None:
            out['number_of_free_fc%d' % order] = int(vs[0])
    
    ### Number of non-zero FCs
    vs = _extract_data(filename, "Number of non-zero  HARMONIC FCs")
    if vs is not None:
        out['number_of_nonzero_fc2'] = int(vs[0])
    #
    for order in range(3, 10):
        vs = _extract_data(filename, "Number of non-zero   ANHARM%d FCs" % order)
        if vs is not None:
            out['number_of_nonzero_fc%d' % order] = int(vs[0])
    
    return out

def read_log_fc2(directory):
    
    filename = directory+'/'+out_dirs['harm']['force']+'/fc2.log'
    return read_log_fc(filename)


#def read_log_fc_anharm(directory):
#
#    out_lasso = read_log_fc_lasso(directory)
#    if out_lasso is None:
#        out_fc3 = read_log_fc3(directory)
#        return out_fc3
#    else:
#        return out_lasso

def read_log_suggest(directory, order=1):
    
    if order == 1:
        mode = 'harm'
    elif order == 2:
        mode = 'cube'    
    
    filename = directory+'/'+out_dirs[mode]['suggest']+'/suggest.log'
    
    ##
    if os.path.exists(filename) == False:
        return None
    
    out = {}
    v = _get_alamode_runtime(filename)
    if v is not None:
        out['time'] = v
    
    nfcs = int(_extract_data(
                filename, "number of  harmonic fcs", back_id=-1
                )[0])
    
    if nfcs is None:
        nfcs = int(_extract_data(
            filename, "number of harmonic fcs", back_id=-1
            )[0])
    
    data = _extract_lines(filename, "space group:")[0].translate(
            str.maketrans({"(": " ", ")": " "})).split()
    
    out['space_group'] = {
            'international': data[-2], 
            'number': int(data[-1]),
            }
    
    out['number_of_symmetry'] = int(_extract_data(
        filename, "number of symmetry operations")[0])
    
    out['number_of_fcs'] = nfcs
    
    out['number_of_structures'] = int(_extract_data(
        filename, "number of disp. patterns")[0])
    
    return out

def read_log_kappa(directory):
    
    outs = {}
    for fc3_type in ['fd', 'lasso']:
        line = directory+'/'+out_dirs['cube']['kappa_%s' % fc3_type] + '_*/kappa.log'
        fns = glob.glob(line)
        for i, fn in enumerate(fns):
            out = read_log_kappa_each(fn)
            label = "%s:%dx%dx%d" % (
                    fc3_type,
                    out['kgrid'][0],
                    out['kgrid'][1],
                    out['kgrid'][2],
                    )
            outs[label] = out
    return outs
   

def read_log_kappa_each(filename):
    
    #filename = directory+'/'+out_dirs['cube']['kappa']+'/kappa.log'
    #
    #if os.path.exists(filename) == False:
    #    filename = directory+'/'+out_dirs['higher']['kappa']+'/kappa.log'
    
    if os.path.exists(filename) == False:
        
        return None
    
    else:
        out = {}
        
        ### volume
        v = _extract_data(
                filename, 
                'volume of the primitive cell', back_id=-2)
        if v is not None:
            out['volume_of_primitive'] = {'value': float(v[0]), 'unit': "A^3"}
        
        out['number_of_atoms_primitive'] = int(
                _extract_data(filename, 'number of atoms in the primitive')[0])
        out['number_of_atoms_supercell'] = int(
                _extract_data(filename, 'number of atoms in the supercell')[0])
        out['number_of_symmetry_operation'] = int(
                _extract_data(filename, 'number of symmetry operations')[0])
        nk1 = int(_extract_data(filename, 'nk1:')[0])
        nk2 = int(_extract_data(filename, 'nk2:')[0])
        nk3 = int(_extract_data(filename, 'nk3:')[0])
        out['kgrid'] = [nk1, nk2, nk3]
        out['number_of_kpoints'] = int(_extract_data(
            filename, "number of k points")[0])
        out['number_of_irreducible_kpoints'] = int(_extract_data(
            filename, "number of irreducible k points")[0])
        v = _get_alamode_runtime(filename)
        if v is not None:
            out['time'] = v
        
        return out

def _get_forces(filename):
    """ Return absolute value of forces 
    Return
    --------
    forces : shape=(natoms,)
    """
    try:
        atoms = ase.io.read(filename, format='vasp-xml')
        all_forces = atoms.get_forces()
        forces = np.zeros(len(atoms))
        for ia in range(len(atoms)):
            forces[ia] = np.linalg.norm(all_forces[ia])
        return forces
    except Exception:
        return None

def _read_each_vaspjob(directory):
    
    file_xml = directory+'/vasprun.xml'
    if os.path.exists(file_xml) == False:
        return None
    
    out = {}
    try:
        out['max_force'] = float(np.max(_get_forces(file_xml)))
    
        ### read outcar
        file_out = directory+'/OUTCAR'
        if os.path.exists(file_out) == False:
            return out
        
        ### time from OUTCAR
        v = _extract_data(file_out, "Total CPU time used")[0]
        if v is not None:
            out['time'] = {'value': float(v), 'unit': 'sec'}
        return out
    
    except Exception:
        return None

def read_log_relax(directory):
    dir_vasp = directory+'/'+out_dirs['relax']
    if os.path.exists(dir_vasp) == False:
        return None
    
    out = {}
    count = 0
    for type in ['full', 'freeze']:
        for i in range(10):
            label = "%s-%d" % (type, i+1)
            
            if type == 'full' and i == 0:
                diri = dir_vasp
            else:
                diri = dir_vasp + "/" + label
            
            fn = diri + '/vasprun.xml'
            if os.path.exists(fn) == False:
                continue
            if count == 0:
                out[type] = _read_each_vaspjob(diri)
            count += 1
        ##
        out['full']['repeat'] = count
    return out

def read_log_nac(directory):
    dir_vasp = directory+'/'+out_dirs['nac']
    out = _read_each_vaspjob(dir_vasp)
    return out

def read_log_forces(directory, mode, fc3_type=None):
    
    ### get the directory name for force calculation
    if mode == 'harm':
        
        dir1 = directory+'/'+out_dirs[mode]['force']
    
    elif mode == 'cube':
        
        dir1 = directory+'/'+out_dirs[mode]['force_%s' % fc3_type]
    
    elif mode == 'higher':
        
        dir1 = directory + '/' + out_dirs[mode]['force']
    
    else:
        warnings.warn(" Warning: %s is not supported." % mode)
        return None

    if os.path.exists(dir1) == False:
        return None
        
    ##
    out = {}
    fmaxes = []
    total_time = 0.
    count = 0
    for i in range(500):
        if i == 0:
            prefix = 'prist'
        else:
            prefix = str(i)
        
        ### read vasprun.xml
        dir_vasp = dir1 + '/'+prefix
        
        if os.path.exists(dir_vasp) == False:
            break
        
        out_each = _read_each_vaspjob(dir_vasp)

        if out_each is None:
            print(" Cannot find %s or the calculation has not been done." % dir_vasp)
            #warnings.warn(" Error in %s" % dir_vasp)
            continue
        
        ## total time
        if 'time' in out_each:
            if out_each['time']['value'] is not None:
                total_time += out_each['time']['value']
        
        ## forces
        if i == 0:
            out['max_force_prist'] = out_each['max_force']
        else:
            count += 1
            fmaxes.append(out_each['max_force'])
       
    out['number_of_patterns'] = count
    out['time'] = {'value': total_time, 'unit': 'sec'}
    
    fmaxes = np.asarray(fmaxes)
    out['minimum_fmax_patterns'] = float(np.min(fmaxes))
    return out

def read_log_lasso(directory):
    
    dir_lasso = directory+'/'+out_dirs['higher']['lasso']
    if os.path.exists(dir_lasso) == False:
        return None
    
    ### force (vasprun.xml and OUTCAR)
    out = {}
    out['force'] = read_log_forces(directory, 'lasso')
    
    ### cv.log
    filename = directory+'/'+out_dirs['higher']['cv']+'/cv.log'
    if os.path.exists(filename):
        out['cv'] = {}
        v = _get_alamode_runtime(filename)
        if v is not None:
            out['cv']['time'] = v
     
    ### lasso.log
    filename = directory+'/'+out_dirs['higher']['lasso']+'/lasso.log' 
    v = read_log_fc(filename)
    if v is not None:
        out['lasso'] = v
    
    return out

#def _analyze_time(out):
#    
#    durations = {}
#    for mode in ['harm', 'cube', 'higher']:
#        if mode in out:
#            if 'force' in out[mode]:
#                if 'time' in out[mode]['force']:
#                    durations['%s_forces' % mode] = \
#                            out[mode]['force']['time']['value']
#    if 'kappa' in out:
#        if 'time' in out['kappa']:
#            durations['kappa'] = out['kappa']['time']['value']

def _get_cellsize_from_log(filename, type=None):
    
    if type.lower() == 'supercell':
        sword = "Supercell"
    elif type.lower() == 'primitive':
        sword = "Primitive"
    
    istart = None
    lines = open(filename, 'r').readlines()
    for il, line in enumerate(lines):
        if sword in line:
            istart = il
            break
    
    if istart is None:
        return None
    ##
    cell = []
    for i in range(3):
        data = lines[istart+2+i].split()
        if len(data) <= 3:
            return None
        else:
            cell.append([float(data[j]) for j in range(3)])
    
    return cell
    
def _get_minimum_frequency(filename):
    
    lines = open(filename, 'r').readlines()
    
    out = {}
    istart = None
    nks = None
    for il, line in enumerate(lines):
        
        if "NONANALYTIC" in line:
            out['nac'] = int(line.split()[2])
        
        if "Number of k points" in line:
            nks = int(line.split()[-1])
            out['number_of_kpoints'] = nks
        
        if "Number of irreducible k points" in line:
            nks = int(line.split()[-1])
            out['number_of_irreducible_kpoints'] = nks
         
        if "Number of atoms in the primitive cell" in line:
            out['number_of_atoms_primitive'] = int(line.split()[-1])
        
        if "Number of atoms in the supercell" in line:
            out['number_of_atoms_supercell'] = int(line.split()[-1])
        
        if "Phonon frequencies below:" in line:
            istart_freq = il
    
    ### get frequencies
    frequencies = []
    out['number_of_bands'] = 3 * out['number_of_atoms_primitive']
    for ik in range(nks):
        for ib in range(out['number_of_bands']):
            num = 4 + istart_freq + ik*(3 + out['number_of_bands']) + ib
            data = lines[num].split()
            frequencies.append(float(data[1]))
    frequencies = np.asarray(frequencies)
    out['minimum_frequency'] = float(np.min(frequencies))
    return out 

def read_log_eigen(directory, mode='band'):
    """ """
    if mode == 'band':
        filename = directory+'/'+out_dirs['harm']['bandos']+'/band.log'
    elif mode == 'evec':
        filename = directory+'/'+out_dirs['harm']['evec']+'/evec_commensurate.log'
    elif mode == 'dos':
        filename = directory+'/'+out_dirs['harm']['bandos']+'/dos.log'
    else:
        warnings.warn(" Error: %s is not supported." % mode)
        return None
    ##
    if os.path.exists(filename) == False:
        return None
    out = {}
    out['supercell'] = _get_cellsize_from_log(filename, type='supercell')
    out['primitive'] = _get_cellsize_from_log(filename, type='primitive')
    v = _get_alamode_runtime(filename)
    if v is not None:
        out['time'] = v
    out.update(_get_minimum_frequency(filename))
    return out

def get_ak_logs(directory):
    """ Return diffent information in log files

    How to Use
    --------------
    >>> 
    >>> directory = './mp-149'
    >>> out_all = get_ak_logs(directory)
    >>> 
    >>> import yaml
    >>> outfile = directory + '/result/log.yaml'
    >>> with open(outfile, "w") as f:
    >>>     yaml.dump(out_all, f)
    >>>     print(" Output", outfile)
    >>> 
    """
    out_all = {"relax":{}, "nac":{}, "harm":{}, "cube":{}, "higher":{}}
    
    ### relax and nac
    v = read_log_relax(directory)
    if v is not None:
        out_all['relax'] = v
    
    v = read_log_nac(directory)
    if v is not None:
        out_all['nac'] = v
    
    ### harmonic
    v = read_log_suggest(directory, order=1)
    if v is not None:
        out_all['harm']["suggest"] = v
    
    v = read_log_forces(directory, 'harm')
    if v is not None:
        out_all['harm']['force'] = v
    
    v = read_log_fc2(directory)
    if v is not None:
        out_all['harm']["fc"] = v
    
    ## harmonic: band, dos, and commensurate points
    for mode in ['band', 'dos', 'evec']:
        v = read_log_eigen(directory, mode=mode)
        if v is not None:
            out_all['harm'][mode] = v
    
    ### cube
    for fc3_type in ['fd', 'lasso']:
        
        v = read_log_forces(directory, 'cube', fc3_type=fc3_type)
        if v is not None:
            out_all['cube']['force_%s' % fc3_type] = v
        
        filename = directory+'/'+out_dirs['cube']['force_%s' % fc3_type]+'/fc3.log' 
        if os.path.exists(filename):
            v = read_log_fc(filename)
            if v is not None:
                out_all['cube']["fc_%s" % fc3_type] = v
    
    ### lasso
    v = read_log_lasso(directory)
    if v is not None:
        out_all['higher'] = v
    
    ### kappa
    try:
        v = read_log_kappa(directory)
        if v is not None:
            out_all["kappa"] = v
    except Exception:
        pass

    ###
    out_mod = {}
    for key in out_all:
        if any(out_all[key]):
            out_mod[key] = out_all[key].copy()
    out_all = out_mod.copy()

    ### total time
    ##times = _analyze_time(out_all)

    return out_all

def get_parser():

    parser = OptionParser()
    
    parser.add_option("-d", "--directory", dest="directory", type="string",
            help="directory name")
    
    parser.add_option("-o", "--outfile", dest="outfile", type="string",
            default=None, help="output .yaml file name")
    
    parser.add_option("-f", "--figname", dest="figname", type="string",
            default=None, help="output figure name")
    
    (options, args) = parser.parse_args()
    
    return options

def main():
    
    options = get_parser()
    
    dirname = options.directory + '/' + out_dirs['result']
    if os.path.exists(dirname) == False:
        print("")
        print(" Cannot find data in %s" % options.directory)
        print("")
        exit()

    log = AkLog(options.directory)
    
    log.write_yaml(outfile=options.outfile) 
    
    log.plot_times(figname=options.figname)

