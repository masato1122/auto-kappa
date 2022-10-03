# -*- coding: utf-8 -*-
import os.path
import numpy as np
from optparse import OptionParser
import ase.io
import datetime
import yaml

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
        
        times = []
        labels = []
        for key in self.get_times():
            
            if '_' in key:
                data = key.split('_')
                lab = "%s(%s)" % (data[1], data[0])
            else:
                lab = key
            
            times.append(self.get_times()[key])
            labels.append(lab)
        
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
                filename, 'job started at')[0].split('at')[1].replace('\n', "")
        line1 = _extract_lines(
                filename, 'job finished at')[-1].split('at')[1].replace("\n", "")
        
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
    else:
        out = {}
        out['time'] = _get_alamode_runtime(filename)
        
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

def read_log_fc3(directory):
    
    filename = directory+'/'+out_dirs['cube']['force']+'/fc3.log' 
    return read_log_fc(filename)

#def read_log_fc_anharm(directory):
#
#    out_lasso = read_log_fc_lasso(directory)
#    if out_lasso is None:
#        out_fc3 = read_log_fc3(directory)
#        return out_fc3
#    else:
#        return out_lasso

def read_log_suggest(directory):
    
    filename = directory+'/'+out_dirs['harm']['suggest']+'/suggest.log'
    if os.path.exists(filename) == False:
        return None
    out = {}
    out['time'] = _get_alamode_runtime(filename)
    nfcs = int(_extract_data(
        filename, "number of  harmonic fcs", back_id=-1)[0])
    if nfcs is None:
        nfcs = int(_extract_data(
            filename, "number of harmonic fcs", back_id=-1)[0])
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
    filename = directory+'/'+out_dirs['cube']['kappa']+'/kappa.log'
    if os.path.exists(filename) == False:
        filename = directory+'/'+out_dirs['lasso']['kappa']+'/kappa.log'
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
        out['nks'] = [nk1, nk2, nk3]
        out['number_of_kpoints'] = int(_extract_data(
            filename, "number of k points")[0])
        out['number_of_irreducible_kpoints'] = int(_extract_data(
            filename, "number of irreducible k points")[0])
        out['time'] = _get_alamode_runtime(filename)
        
        return out

def _get_forces(filename):
    """ Return absolute value of forces 
    Return
    --------
    forces : shape=(natoms,)
    """
    atoms = ase.io.read(filename, format='vasp-xml')
    all_forces = atoms.get_forces()
    forces = np.zeros(len(atoms))
    for ia in range(len(atoms)):
        forces[ia] = np.linalg.norm(all_forces[ia])
    return forces

def _read_each_vaspjob(directory):
    
    file_xml = directory+'/vasprun.xml'
    if os.path.exists(file_xml) == False:
        return None
    
    out = {}
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

def read_log_relax(directory):
    dir_vasp = directory+'/'+out_dirs['relax']
    if os.path.exists(dir_vasp) == False:
        return None
    out = _read_each_vaspjob(dir_vasp)
    return out

def read_log_nac(directory):
    dir_vasp = directory+'/'+out_dirs['nac']
    out = _read_each_vaspjob(dir_vasp)
    return out

def read_log_forces(directory, mode):
    
    dir1 = directory+'/'+out_dirs[mode]['force']
    if os.path.exists(dir1) == False:
        return None
    
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
    
    dir_lasso = directory+'/'+out_dirs['lasso']['lasso']
    if os.path.exists(dir_lasso) == False:
        return None
    
    ### force (vasprun.xml and OUTCAR)
    out = {}
    out['force'] = read_log_forces(directory, 'lasso')
    
    ### cv.log
    filename = directory+'/'+out_dirs['lasso']['cv']+'/cv.log'
    if os.path.exists(filename):
        out['cv'] = {}
        time = _get_alamode_runtime(filename)
        if time is not None:
            out['cv']['time'] = time
    
    ### lasso.log
    filename = directory+'/'+out_dirs['lasso']['lasso']+'/lasso.log' 
    v = read_log_fc(filename)
    if v is not None:
        out['lasso'] = v
    
    return out

def _analyze_time(out):
    
    durations = {}
    for mode in ['harm', 'cube', 'lasso']:
        if mode in out:
            if 'force' in out[mode]:
                if 'time' in out[mode]['force']:
                    durations['%s_forces' % mode] = \
                            out[mode]['force']['time']['value']
    if 'kappa' in out:
        if 'time' in out['kappa']:
            durations['kappa'] = out['kappa']['time']['value']

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
    out_all = {"relax":{}, "nac":{}, "harm":{}, "cube":{}, "lasso":{}}
     
    ### relax and nac
    v = read_log_relax(directory)
    if v is not None:
        out_all['relax'] = v
    
    v = read_log_nac(directory)
    if v is not None:
        out_all['nac'] = v
    
    ### harmonic
    v = read_log_suggest(directory)
    if v is not None:
        out_all['harm']["suggest"] = v
    
    v = read_log_forces(directory, 'harm')
    if v is not None:
        out_all['harm']['force'] = v
    
    v = read_log_fc2(directory)
    if v is not None:
        out_all['harm']["fc"] = v
    
    ### cube
    v = read_log_forces(directory, 'cube')
    if v is not None:
        out_all['cube']['force'] = v
    
    v = read_log_fc3(directory)
    if v is not None:
        out_all['cube']["fc"] = v
    
    ### lasso
    v = read_log_lasso(directory)
    if v is not None:
        out_all['lasso'] = v
    
    ### kappa
    v = read_log_kappa(directory)
    if v is not None:
        out_all["kappa"] = v
    
    ###
    out_mod = {}
    for key in out_all:
        if any(out_all[key]):
            out_mod[key] = out_all[key].copy()
    out_all = out_mod.copy()

    ### total time
    times = _analyze_time(out_all)

    return out_all

def main(options):
    
    out_all = get_ak_logs(options.directory)
    
    with open(options.outfile, "w") as f:
        yaml.dump(out_all, f)
        print(" Output", options.outfile)


if __name__ == '__main__':
    parser = OptionParser()
    
    parser.add_option("-d", "--directory", dest="directory", type="string",
            help="input file name")
    
    parser.add_option("-o", "--outfile", dest="outfile", type="string",
            default="log.yaml", help="output .yaml file name")
    
    (options, args) = parser.parse_args()
    main(options)

