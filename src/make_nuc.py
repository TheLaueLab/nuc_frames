from h5py import File, Group
import numpy as np
import os, random, time

def import_contacts(filePath):
  
  fileObj = open(filePath, 'r')
  
  contactDict = {}
   
  for line in fileObj:
    chrA, f_start_a, f_end_a, start_a, end_a, strand_a, chrB, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_group, pair_id, swap_pair = line.split()
    
    if strand_a == '+':
      posA = int(f_start_a)
    else:
      posA = int(f_end_a)
    
    if strand_b == '+':
      posB = int(f_start_b)       
    else:
      posB = int(f_end_b)
 
    if chrA > chrB:
      chrA, chrB = chrB, chrA
      posA, posB = posB, posA
 
    chromoPair = (chrA, chrB)
 
    if chromoPair in contactDict:
      contactDict[chromoPair].append((posA, posB, 1, ambig_group))
    else:
      contactDict[chromoPair] = [(posA, posB, 1, ambig_group)]
   
  fileObj.close()
  
  for chromoPair in contactDict:
    contactDict[chromoPair] = np.array(contactDict[chromoPair], np.uint32)
  
  return contactDict


def import_coords(filePath):
    
  fileObj = open(filePath, 'r')
  
  posDict = {}
  coordsDict = {}  
  chromo = None
  
  for line in fileObj:
    
    data = line.split()
    nItems = len(data)
    
    if not nItems:
      continue
    
    elif data[0] == '#':
      continue
    
    elif nItems == 3:
      chromo, nCoords, nModels = data
      nCoords = int(nCoords)
      nModels = int(nModels)
      
      chromoPos = []
      chromoCoords = np.empty((nModels, nCoords, 3), np.float32)
      
      coordsDict[chromo] = chromoCoords
      posDict[chromo] = chromoPos
      
      nCoords = int(nCoords)
      nModels = int(nModels)
      check = (nModels * 3) + 1
      i = 0
      
    elif not chromo:
      raise Exception('Missing chromosome record in file %s' % filePath)
     
    elif nItems != check:
      msg = 'Data size in file %s does not match Position + Models * Positions * 3'
      raise Exception(msg % filePath)
    
    else:
      chromoPos.append(int(data[0]))
      
      coord = [float(x) for x in data[1:]]
      coord = np.array(coord).reshape(nModels, 3)
      chromoCoords[:,i] = coord
      i += 1
  
  fileObj.close()
  
  return posDict, coordsDict



def make_nuc(ncc_file_path, n3d_file_path, out_file_name):
  
  if not out_file_name.lower().endswith('.nuc'):
    out_file_name = out_file_name + '.nuc'
  
  contact_dict = import_contacts(ncc_file_path)
  
  contact_name = os.path.splitext(os.path.basename(ncc_file_path))[0]
  
  pos_dict, coords_dict = import_coords(n3d_file_path)
  
  root = File(out_file_name, mode='w')
        
  hierarchy = (('contacts',   ('original', 'working')),
                ('display',    ()),
                ('chromosomes',()),
                ('dataTracks', ('derived', 'external', 'innate')),
                ('sample',     ('protocol', 'organism', 'tissue')),
                ('structures', ('0')),
                ('images',     ())
                )
   
  for parent, children in hierarchy:
    group = root.create_group(parent)
  
    for child in children:
      group.create_group(child)
  
  for child in ('particles', 'restraints', 'transforms', 'coords'):
    root['structures']['0'].create_group(child)
  
  now = int(time.time())
  random.seed(now)        
  
  root.attrs['id'] = np.array([random.random(), now, now], np.float32)
  
  root['sample'].attrs['name'] = np.string_('Unknown')  
  
  contact_group = root['contacts']['working'].create_group(contact_name)
  
  for chromoPair in contact_dict:
    chrA, chrB = chromoPair
    
    if chrA not in contact_group:
      contact_group.create_group(chrA)

    contact_group[chrA].create_dataset(chrB, dtype=np.uint32, data=contact_dict[chromoPair].T)
    
  coords_group   = root['structures']['0']['coords']
  particle_group = root['structures']['0']['particles']
 
  
  for chromo in coords_dict:
    coords_group.create_dataset(chromo, dtype=np.float64, data=coords_dict[chromo])
    
    pos = np.array(pos_dict[chromo], np.uint32)
    group = particle_group.create_group(chromo)
    group.create_dataset('positions', dtype=np.uint32, data=pos)
    
    chromo_group = root['chromosomes'].create_group(chromo)
    chromo_group.attrs['limits'] = np.array([pos.min(), pos.max()])
    
  root.flush()
  
  
  
if __name__ == '__main__':

  from argparse import ArgumentParser
    
  epilog = 'For further help on running this program please email tjs23@cam.ac.uk or wb104@cam.ac.uk'
  
  arg_parse = ArgumentParser(prog='make_nuc', description='Converter to make legacy .nuc files from .ncc (contact) and .n3d (3D coorrinate) files',
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('ncc_path', nargs=1, metavar='NCC_FILE',
                         help='Input NCC format file containing single-cell Hi-C contact data')

  arg_parse.add_argument('n3d_path', nargs=1, metavar='N3D_FILE',
                         help='Input N3D format file containing genome structure coordinates')

  arg_parse.add_argument('nuc_path', nargs=1, metavar='NUC_FILE',
                         help='Output file path for .nuc format file')
  
  args = vars(arg_parse.parse_args())
  
  ncc_file_path = args['ncc_path'][0]
  
  n3d_file_path = args['n3d_path'][0]
  
  nuc_file_path = args['nuc_path'][0]
  
  make_nuc(ncc_file_path, n3d_file_path, nuc_file_path)
  
  
  

    
  
