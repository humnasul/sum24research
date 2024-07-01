import h5py

def print_hdf5_data(h5_file):
  """Prints the data in an HDF5 file.

  Args:
    h5_file: The HDF5 file to print.
  """

  for group in h5_file.keys():
    print(group)
    for dataset in h5_file[group].keys():
      print(dataset)
      print(h5_file[group][dataset])

if __name__ == '__main__':
  with h5py.File('A20133370505_BLD4_L2.hdf', 'r') as f:
    print_hdf5_data(f)