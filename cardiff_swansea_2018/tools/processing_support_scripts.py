from thetis import *
import pickle

def output_field_h5(output_directory, field, name):
    checkpoint_file = checkpointing.DumbCheckpoint(output_directory + "/" + name)
    checkpoint_file.store(field)
    checkpoint_file.close()