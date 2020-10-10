###################################################################################################
#    This file is part of parKVFinder.                                                            #
#                                                                                                 #
#    parKVFinder is free software: you can redistribute it and/or modify                          #
#    it under the terms of the GNU General Public License as published by                         #
#    the Free Software Foundation, either version 3 of the License, or                            #
#    (at your option) any later version.                                                          #
#                                                                                                 #
#    parKVFinder is distributed in the hope that it will be useful,                               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                               #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                #
#    GNU General Public License for more details.                                                 #
#                                                                                                 #
#    You should have received a copy of the GNU General Public License                            #
#    along with parKVFinder.  If not, see <http://www.gnu.org/licenses/>.                         #
#                                                                                                 #
###################################################################################################

###################################################################################################
#    This is the source code of the PyMOL parKVFinder Tools Interface for PyMOL. It was           #
#    developed using tk and Python.                                                               #
#                                                                                                 #
#    After the installation of parKVFinder, a copy of this file will be in the tools directory.   #
#    Changes in this file are not advised, as it controls all of the parKVFinder features         #
#    in PyMOL.                                                                                    #
###################################################################################################

# Import required modules
import os
import subprocess
import math
import time
import Pmw
from pymol import cmd
import pymol
import toml
import sys

try:
    from builtins import chr  # unichr works for python 2 and 3
except ImportError:
    from __builtin__ import chr

# Import Tkinter
try:
    import Tkinter  # Python 2
    from Tkinter import *  # Python 2
    import tkFileDialog
    import tkMessageBox
except ImportError:
    import tkinter as Tkinter  # Python 3
    from tkinter import *  # Python 3
    import tkinter.filedialog as tkFileDialog  # Python 3
    import tkinter.messagebox as tkMessageBox  # Python 3


def __init_plugin__(self):
    self.menuBar.addmenuitem('Plugin',
                             'command',
                             'PyMOL parKVFinder Tools v1.0',
                             label='PyMOL parKVFinder Tools v1.0',
                             command=lambda s=self: PyMOL_parKVFinder_Tools(s))


class PyMOL_parKVFinder_Tools:
    """
    Graphical User Interface
    """
    def __init__(self, app):
        self.parent = app.root

        ################################################################################
        #                                                                              #
        #    Attributes                                                                #
        #                                                                              #
        ################################################################################

        # Get KVFinder_PATH
        KVFinder_PATH = os.getenv('KVFinder_PATH')

        #### MessageBox Warnings #######################################################################################
        # KVFinder_PATH not found by os.getenv()

        #### Check KVFinder_PATH in .bash_profile
        if KVFinder_PATH is None:
            if os.path.exists(os.getenv('HOME')+"/.bash_profile"):
                file_env = open(os.getenv('HOME')+"/.bash_profile")
                for line in file_env:
                    if line.find('export KVFinder_PATH') == 0:
                        KVFinder_PATH = line.split("=")[1].rstrip('\n')
            # MessageBox Warning 1: Found KVFinder_PATH in .bash_profile
            if KVFinder_PATH is not None:
                tkMessageBox.showwarning("Warning",
                                         """Check File Locations!\n \
KVFinder_PATH was found in .bash_profile file.""",
                                         parent=self.parent)

        #### Check KVFinder_PATH in .bashrc
        if KVFinder_PATH is None:
            if os.path.exists(os.getenv('HOME') + "/.bashrc"):
                file_env = open(os.getenv('HOME') + "/.bashrc")
                for line in file_env:
                    if line.find('export KVFinder_PATH') == 0:
                        KVFinder_PATH = line.split("=")[1].rstrip('\n')
            if KVFinder_PATH is not None:
                tkMessageBox.showwarning("Warning",
                                         """Check File Locations!\n KVFinder_PATH was found in .bashrc file.""",
                                         parent=self.parent)

        #### KVFinder_PATH were not found
        if KVFinder_PATH is None:
            tkMessageBox.showwarning("Warning",
                                     """KVFinder_PATH not found: set File Locations! \n\
Otherwise, parKVFinder cannot be executed in PyMOL parKVFinder Tools.""",
                                     parent=self.parent)
            # Set an empty KVFinder_PATH
            KVFinder_PATH = ""

        # Default values for attributes
        self.defaults = {
            # Environment variables
            'KVFinder_PATH': "" if KVFinder_PATH is None else KVFinder_PATH,
            # Program executables
            'parKVFinder': "" if KVFinder_PATH is None else "{}/{}".format(str(KVFinder_PATH), 'parKVFinder'),
            # Parameters
            'probe_in': 1.4,
            'probe_out': 4.0,
            'step_size': 0.0,
            'resolution': 'Low',
            'volume_cutoff': 5.0,
            'ligand_cutoff': 5.0,
            'removal_distance': 2.4,
            'base_name': 'output',
            # Modes
            'grid_method': False,
            'ligand_mode': False,
            'search_procedure': 'Whole Protein',
            'surface_representation': "Molecular Surface (VdW)",
            'cavity_representation': "Filtered (less memory consumption)",
            # Files
            'dictionary': "" if KVFinder_PATH is None else "{}/{}".format(str(KVFinder_PATH), 'dictionary'),
        }

        # Get working directory
        self.wd = StringVar()
        self.wd.set(os.getcwd())

        # Set modes and methods to default
        self.grid_method = BooleanVar()
        self.grid_method.set(self.defaults['grid_method'])
        self.ligand_mode = BooleanVar()
        self.ligand_mode.set(self.defaults['ligand_mode'])
        self.loaded_results = BooleanVar()
        self.loaded_results.set(0)

        # Declare results variables
        # Set last output object
        self.last_object = StringVar()
        self.last_object.set("")
        # Set last input pdb
        self.last_input = StringVar()
        self.last_input.set("")
        # Set last ligand pdb
        self.last_ligand = StringVar()
        self.last_ligand.set("")

        # Print working directory
        sys.stdout.write('Working directory: {}/'.format(self.wd.get()))

        ################################################################################
        #                                                                              #
        #    (0) Dialog Interface                                                      #
        #                                                                              #
        ################################################################################

        # Create the dialog.
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=('Run parKVFinder',
                                          'Show Grid',
                                          'Save Parameters',
                                          'Restore Default Values',
                                          'Exit'),
                                 defaultbutton='Run parKVFinder',
                                 title='PyMOL parKVFinder Tools v1.0',
                                 command=self.execute)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        # Create Dialog label
        self.w = Label(self.dialog.interior(),
                       text=
                       """\nparKVFinder (parallel KVFinder) software identifies and describes cavities in a target \
biomolecular structure using a dual probe system.\n\nThe description includes spatial and constitutional \
characterization. The spatial description includes shape, volume and area. The constitutional description includes \
amino acids that form the identified cavities.\n""",
                       background='lightgray',
                       foreground='black',
                       justify=LEFT,
                       wraplength=1000,
                       padx=10
                       )
        self.w.pack(side=TOP, fill=X, padx=10, pady=10)

        # Create notebook tabs
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill=BOTH, expand=1, padx=10, pady=5)

        ################################################################################################################
        #    (1) Main Tab                                                                                              #
        ################################################################################################################

        # Set up Main tab
        main = self.notebook.add('Main')

        #### (1.1) Main Options ########################################################################################
        group = Pmw.Group(main,
                          tag_text='Main Options')
        group.pack(fill=BOTH, expand=1, padx=10, pady=5)

        # Create two frames for Main Options
        main_options_frame1 = Frame(group.interior())
        main_options_frame2 = Frame(group.interior())
        main_options_frame1.grid(row=0, column=0)
        main_options_frame2.grid(row=0, column=1)

        #### (1.1.1) Input PDB Group
        input_pdb_group = Pmw.Group(main_options_frame1,
                                    tag_text='Input PDB')
        input_pdb_group.grid(row=0, column=0, padx=20, pady=10)

        ## (1.1.1.1) Input PDB Frame
        input_pdb_frame = Frame(input_pdb_group.interior())
        input_pdb_frame.grid(row=0, padx=10, pady=10)
        scrollbar_input_pdb = Scrollbar(input_pdb_frame,
                                        orient='vertical')
        self.pdb_file = Listbox(input_pdb_frame,
                                height=10,
                                width=40,
                                yscrollcommand=scrollbar_input_pdb.set,
                                exportselection=0)
        scrollbar_input_pdb.config(command=self.pdb_file.yview)
        scrollbar_input_pdb.pack(side=RIGHT, fill=Y)
        self.pdb_file.pack(side=LEFT, fill=BOTH)

        ## (1.1.1.2) Buttons Input PDB
        input_pdb_buttons = Frame(input_pdb_group.interior())
        input_pdb_buttons.grid(row=1)
        input_pdb_refresh_button = Button(input_pdb_buttons,
                                          text="Refresh List",
                                          command=lambda: self.refresh_listbox(self.pdb_file))
        input_pdb_refresh_button.pack(side=LEFT, anchor='n')
        input_pdb_upload_button = Button(input_pdb_buttons,
                                         text="Upload PDB",
                                         command=lambda: self.upload_pdb(self.pdb_file))
        input_pdb_upload_button.pack(side=LEFT, anchor='n')

        #### (1.1.2) Surface Representation Group
        surface_representation_group = Pmw.Group(main_options_frame1,
                                                 tag_text='Surface Representation')
        surface_representation_group.grid(row=1, column=0, padx=20, pady=10)

        ## (1.1.2.1) Surface Representation Radio Button
        self.surface_representation = {"Solvent Accessible Surface (SAS)": False,
                                       "Molecular Surface (VdW)": True}
        self.surface = Pmw.RadioSelect(surface_representation_group.interior(),
                                       buttontype='radiobutton',
                                       orient='vertical',
                                       labelpos='w')
        for text in self.surface_representation.keys():
            self.surface.add(text)
        self.surface.setvalue(self.defaults['surface_representation'])
        self.surface.pack(padx=10, pady=5, fill=X)

        #### (1.1.3) Settings Group
        settings_group = Pmw.Group(main_options_frame2,
                                   tag_text='Settings')
        settings_group.grid(row=0, padx=20, pady=10)

        ## (1.1.3.1) Base Name
        self.base_name = Pmw.EntryField(settings_group.interior(),
                                        labelpos='w',
                                        label_text='Output Base Name:',
                                        value=self.defaults['base_name'],
                                        entry_width=6)
        self.base_name.pack(anchor='w', fill=X, padx=10)

        ## (1.1.3.2) Probe In
        self.probe_in = Pmw.Counter(settings_group.interior(),
                                    labelpos='w',
                                    label_text=u"Probe In ({}):".format(chr(0x212b)),
                                    entryfield_value=self.defaults['probe_in'],
                                    entryfield_validate={'validator': 'real', 'min': 0.0, 'max': 5.0},
                                    datatype={'counter': self._custom_real_counter},
                                    entry_width=6,
                                    increment=0.10)
        self.probe_in.pack(anchor='w', fill=X, padx=10)

        ## (1.1.3.3) Probe Out
        self.probe_out = Pmw.Counter(settings_group.interior(),
                                     labelpos='w',
                                     label_text=u"Probe Out ({}):".format(chr(0x212b)),
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': 0.0, 'max': 50.0},
                                     increment=0.10,
                                     entry_width=6,
                                     entryfield_value=self.defaults['probe_out'])
        self.probe_out.pack(anchor='w', fill=X, padx=10)

        ## (1.1.3.4) Volume Cutoff
        self.volume_cutoff = Pmw.Counter(settings_group.interior(),
                                         labelpos='w',
                                         label_text=u"Volume Cutoff ({}{}):".format(chr(0x212b), chr(0x00b3)),
                                         datatype={'counter': self._custom_real_counter},
                                         entryfield_validate={'validator': 'real', 'min': 0.0},
                                         increment=0.10,
                                         entry_width=6,
                                         entryfield_value=self.defaults['volume_cutoff'])
        self.volume_cutoff.pack(anchor='w', fill=X, padx=10)

        ## (1.1.3.5) Removal Distance
        self.removal_distance = Pmw.Counter(settings_group.interior(),
                                            labelpos='w',
                                            label_text=u"Removal Distance ({})".format(chr(0x212b)),
                                            datatype={'counter': self._custom_real_counter},
                                            entryfield_validate={'validator': 'real', 'min': 0.0, 'max': 10.0},
                                            increment=0.10,
                                            entry_width=6,
                                            entryfield_value=self.defaults['removal_distance'])
        self.removal_distance.pack(anchor='w', fill=X, padx=10)

        # Align labels
        bars = (self.probe_in, self.probe_out, self.volume_cutoff, self.removal_distance)
        Pmw.alignlabels(bars)

        ## (1.1.3.6) Grid Spacing Methods Group
        grid_settings_group = Pmw.Group(settings_group.interior(),
                                        tag_text='Grid Spacing Methods')
        grid_settings_group.pack(fill=X, pady=10, padx=10)

        # Create Grid Method Options
        grid_groups = []
        radio_grid_settings_frame = Frame(grid_settings_group.interior())
        radio_grid_settings_frame.pack(padx=10, pady=10, expand=1, fill=BOTH)

        # (1.1.3.6.1) Step Size (Explicit method)
        grid_methods = Pmw.Group(radio_grid_settings_frame,
                                 tag_pyclass=Radiobutton,
                                 tag_text=u"Step Size ({}):".format(chr(0x212b)),
                                 tag_value=1,
                                 tag_variable=self.grid_method,
                                 tag_command=self.check_grid_method)
        grid_methods.pack(fill=BOTH, expand=1, side=LEFT)
        grid_groups.append(grid_methods)
        self.step_size = Pmw.Counter(grid_methods.interior(),
                                     labelpos='w',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': 0.0, 'max': 2.0},
                                     increment=0.1,
                                     entry_width=10,
                                     entryfield_value=self.defaults['step_size'])
        self.step_size.pack(anchor='w')

        # (1.1.3.6.2) Resolution (Implicit method)
        grid_methods = Pmw.Group(radio_grid_settings_frame,
                                 tag_pyclass=Radiobutton,
                                 tag_text="Resolution",
                                 tag_value=0,
                                 tag_variable=self.grid_method,
                                 tag_command=self.check_grid_method)
        grid_methods.pack(fill=BOTH, expand=1, side=LEFT)
        self.resolution = Pmw.OptionMenu(grid_methods.interior(),
                                         labelpos='w',
                                         items=('High',
                                                'Medium',
                                                'Low',
                                                'Off'),
                                         command=self.check_step_size_status,
                                         menubutton_width=10,
                                         initialitem=self.defaults['resolution'])
        self.resolution.pack(anchor='w')
        grid_groups.append(grid_methods)

        # Align labels
        Pmw.aligngrouptags(grid_groups)

        #### (1.1.4) Cavities Representation Group
        cavity_representation_group = Pmw.Group(main_options_frame2,
                                                tag_text='Cavity Representation')
        cavity_representation_group.grid(row=1, padx=0, pady=0)

        ## (1.1.4.1) Cavity Representation Radio Button
        self.cavity_representation_options = {"Full": True,
                                              "Filtered (less memory consumption)": False}
        self.cavity_representation = Pmw.RadioSelect(cavity_representation_group.interior(),
                                                     buttontype='radiobutton',
                                                     orient='vertical',
                                                     labelpos='w')
        for text in sorted(self.cavity_representation_options.keys(), reverse=True):
            self.cavity_representation.add(text)
        self.cavity_representation.setvalue(self.defaults['cavity_representation'])
        self.cavity_representation.pack(padx=10, pady=5, fill=X)

        #### (1.2) File Locations Group ################################################################################
        file_locations_group = Pmw.Group(main, tag_text='File Locations')
        file_locations_group.pack(fill=BOTH, expand=1, padx=10, pady=5)

        #### (1.2.1) parKVFinder Executable Path
        file_locations_group_row1 = Frame(file_locations_group.interior())
        file_locations_group_row1.grid(row=0, column=0, padx=20, pady=5, sticky='w')

        self.execKVFinder = Entry(file_locations_group_row1, takefocus=0, width=80)
        Button(file_locations_group_row1,
               width=15,
               text="parKVFinder:",
               command=lambda: self.select_file(self.execKVFinder,
                                                'Select parKVFinder executable')).pack(side=LEFT)
        self.execKVFinder.insert(0, self.defaults['parKVFinder'])
        self.execKVFinder.configure(state='readonly')
        self.execKVFinder.pack(fill=X, side=LEFT)

        #### (1.2.2) Van der Waals Dictionary Path
        file_locations_group_row2 = Frame(file_locations_group.interior())
        file_locations_group_row2.grid(row=1, column=0, columnspan=3, padx=20, pady=5, sticky='w')

        self.dictionary = Entry(file_locations_group_row2, takefocus=0, width=80)
        Button(file_locations_group_row2,
               width=15,
               text="Dictionary:",
               command=lambda: self.select_file(self.dictionary,
                                                'Select radii dictionary')).pack(side=LEFT)
        self.dictionary.insert(0, self.defaults['dictionary'])
        self.dictionary.configure(state='readonly')
        self.dictionary.pack(fill=X, side=LEFT)

        #### (1.2.3) Output Directory Path
        file_locations_group_row3 = Frame(file_locations_group.interior())
        file_locations_group_row3.grid(row=2, column=0, columnspan=3, padx=20, pady=5, sticky='w')

        self.output = Entry(file_locations_group_row3, takefocus=0, width=80)
        Button(file_locations_group_row3,
               width=15,
               text="Output Directory:",
               command=lambda: self.select_directory(self.output,
                                                     'Select output directory')).pack(side=LEFT)
        self.output.insert(0, self.wd.get())
        self.output.configure(state='readonly')
        self.output.pack(fill=X)

        ################################################################################################################
        #    (2) Search Space Tab                                                                                      #
        ################################################################################################################

        # Set up Search Space tab
        search_space = self.notebook.add('Search Space')

        #### (2.1) Adjustment Methods Group ############################################################################
        adjustment_method_group = Pmw.Group(search_space,
                                            tag_text='Adjustment Methods')
        adjustment_method_group.pack(fill=BOTH, expand=1, padx=10, pady=5)

        # Create frames in Adjustment Methods Group
        self.adjustment_method_group_frame1 = Frame(adjustment_method_group.interior())
        self.adjustment_method_group_frame1.pack(fill=BOTH, side=LEFT)
        adjustment_method_group_frame2 = Frame(adjustment_method_group.interior())
        adjustment_method_group_frame2.pack(fill=BOTH, side=LEFT)

        #### (2.1.1) Search Procedure Group
        search_procedure_group = Pmw.Group(self.adjustment_method_group_frame1,
                                           tag_text='Search Procedure')
        search_procedure_group.pack(padx=10, pady=5, fill=X)

        ## (2.1.1.1) Search Procedure Radio Button
        search_procedure_group_frame = Frame(search_procedure_group.interior(), width=350)
        search_procedure_group_frame.pack(fill=X)
        self.search_procedure_options = {'Whole Protein': True,
                                         'Box Adjustment': False}
        self.search_procedure = Pmw.RadioSelect(search_procedure_group_frame,
                                                buttontype='radiobutton',
                                                orient='vertical',
                                                labelpos='w',
                                                command=self.update_box_adjustment_frame)
        for text in self.search_procedure_options.keys():
            self.search_procedure.add(text)
        self.search_procedure.setvalue('Whole Protein')
        self.search_procedure.pack(padx=10, pady=5, fill=X)

        #### (2.1.2) Box Adjustment Dummy Frame
        self.box_adjustment_frame = Frame(self.adjustment_method_group_frame1, width=350)
        self.box_adjustment_frame.pack(padx=10, pady=5, fill=X)

        #### (2.1.3) Ligand Adjustment Group
        ligand_adjustment_group = Pmw.Group(adjustment_method_group_frame2,
                                            tag_pyclass=Checkbutton,
                                            tag_text='Ligand Adjustment',
                                            tag_variable=self.ligand_mode,
                                            tag_command=self.check_ligand_mode)
        ligand_adjustment_group.pack(fill=BOTH, padx=10, pady=5)

        ## (2.1.3.1) Ligand PDB Group
        ligand_pdb_group = Pmw.Group(ligand_adjustment_group.interior(),
                                     tag_text='Ligand PDB')
        ligand_pdb_group.grid(row=0, column=0, padx=20, pady=10)

        # (2.1.3.1.1) Ligand PDB Frame
        ligand_pdb_group_frame = Frame(ligand_pdb_group.interior())
        ligand_pdb_group_frame.grid(row=0, padx=10, pady=10)
        scrollbar_ligand_pdb = Scrollbar(ligand_pdb_group_frame,
                                         orient="vertical")
        self.ligand_file = Listbox(ligand_pdb_group_frame,
                                   height=10,
                                   width=40,
                                   yscrollcommand=scrollbar_ligand_pdb.set,
                                   exportselection=0)
        scrollbar_ligand_pdb.config(command=self.ligand_file.yview)
        scrollbar_ligand_pdb.pack(side=RIGHT, fill=Y)
        self.ligand_file.pack(side=LEFT, fill=BOTH)

        # (2.1.3.1.2) Button Ligand PDB
        ligand_pdb_buttons = Frame(ligand_pdb_group.interior())
        ligand_pdb_buttons.grid(row=1)

        # Refresh List Button
        self.ligand_pdb_refresh_button = Button(ligand_pdb_buttons,
                                                text="Refresh List",
                                                command=lambda: self.refresh_listbox(self.ligand_file))
        self.ligand_pdb_refresh_button.pack(side=LEFT, anchor='n')
        # Upload PDB Button
        self.ligand_pdb_upload_button = Button(ligand_pdb_buttons,
                                               text="Upload Ligand",
                                               command=lambda: self.upload_pdb(self.ligand_file))
        self.ligand_pdb_upload_button.pack(side=LEFT, anchor='n')
        # Disable Buttons
        self.ligand_pdb_refresh_button.configure(state='disabled')
        self.ligand_pdb_upload_button.configure(state='disabled')

        ## (2.1.3.2) Ligand Cutoff
        self.ligand_cutoff = Pmw.Counter(ligand_adjustment_group.interior(),
                                         labelpos='w',
                                         label_text=u"Ligand Cutoff ({}):".format(chr(0x212b)),
                                         datatype={'counter': self._custom_real_counter},
                                         entryfield_validate={'validator': 'real', 'min': 0.0},
                                         increment=0.10,
                                         entry_width=20,
                                         entryfield_value=self.defaults['ligand_cutoff'])
        self.ligand_cutoff.grid(row=1, pady=10)

        ################################################################################################################
        #    (3) Results Tab                                                                                           #
        ################################################################################################################

        # Set up Results Tab
        results = self.notebook.add('Results Visualization')

        #### (3.1) Information Group ###################################################################################
        information_group = Pmw.Group(results,
                                      tag_text='Information')
        information_group.pack(fill=BOTH, padx=10, pady=5)

        #### (3.1.1) Results File Frame
        results_frame = Frame(information_group.interior())
        results_frame.grid(row=0, sticky='W', padx=10, pady=1)

        ## (3.1.1.1) Results File Button
        results_file_button = Button(results_frame,
                                     width=8,
                                     text="Results File:",
                                     command=lambda: self.select_file(self.results_file,
                                                                      'Select KVFinder Results File',
                                                                      [
                                                                          ("KVFinder Results File",
                                                                           "*KVFinder.results.toml"
                                                                           ),
                                                                          ('All files',
                                                                           '*')
                                                                      ]
                                                                      )
                                     )
        results_file_button.pack(side=LEFT)

        ## (3.1.1.2) Results File Path
        self.results_file = Entry(results_frame,
                                  takefocus=0,
                                  width=100)
        self.results_file.configure(state='readonly')
        self.results_file.pack(fill=X, side=LEFT)

        ## (3.1.1.3) Load Results File Button
        load_results_button = Button(results_frame,
                                     width=4,
                                     text="Load",
                                     command=lambda: self.load_results(1))
        load_results_button.pack(side=LEFT)

        #### (3.1.2) Files Paths and Step Size
        ## (3.1.2.1) Input File Path
        self.inp_file = Label(information_group.interior(),
                              text="Input File: ")
        self.inp_file.grid(row=1, sticky='W', padx=10, pady=1)
        ## (3.1.2.2) Ligand File Path
        self.lig_file = Label(information_group.interior(),
                              text="Ligand File: ")
        self.lig_file.grid(row=2, sticky='W', padx=10, pady=1)
        ## (3.1.2.3) Output File Path
        self.output_file = Label(information_group.interior(),
                              text="Output File: ")
        self.output_file.grid(row=3, sticky='W', padx=10, pady=1)
        ## (3.1.2.4) Step Size
        self.step_size_info = Label(information_group.interior(),
                                    text=u"Step Size ({}):".format(chr(0x212b)))
        self.step_size_info.grid(row=4, sticky='W', padx=10, pady=1)

        #### (3.2) Descriptors Group ###################################################################################
        descriptors = Pmw.Group(results,
                                tag_text='Descriptors')
        descriptors.pack(fill=BOTH, expand=1, padx=10, pady=5)

        #### (3.2.1) Volume Listbox
        self.volume_frame = Frame(descriptors.interior())
        self.volume_frame.pack(side=LEFT, padx=10)
        self.scrollbar_volume_vertical = Scrollbar(self.volume_frame,
                                                   orient="vertical")
        self.scrollbar_volume_horizontal = Scrollbar(self.volume_frame,
                                                     orient="horizontal")
        Label(self.volume_frame,
              text=u"Volume ({}{})".format(chr(0x212b),
                                           chr(0x00b3))).pack()
        self.list_volume = Listbox(self.volume_frame,
                                   height=18,
                                   width=25,
                                   yscrollcommand=self.scrollbar_volume_vertical.set,
                                   xscrollcommand=self.scrollbar_volume_horizontal.set,
                                   selectmode='multiple',
                                   exportselection=0)
        self.scrollbar_volume_vertical.config(command=self.list_volume.yview)
        self.scrollbar_volume_vertical.pack(side=RIGHT, fill=Y)
        self.scrollbar_volume_horizontal.config(command=self.list_volume.xview)
        self.scrollbar_volume_horizontal.pack(side=BOTTOM, fill=X)
        self.list_volume.pack(side=LEFT, fill=BOTH, expand=1)

        #### (3.2.2) Area Listbox
        self.area_frame = Frame(descriptors.interior())
        self.area_frame.pack(side=LEFT, padx=10)
        self.scrollbar_area_vertical = Scrollbar(self.area_frame,
                                                 orient="vertical")
        self.scrollbar_area_horizontal = Scrollbar(self.area_frame,
                                                   orient="horizontal")
        Label(self.area_frame,
              text=u"Surface Area ({}{})".format(chr(0x212b),
                                                 chr(0x00b2))).pack()
        self.list_area = Listbox(self.area_frame,
                                 height=18,
                                 width=25,
                                 yscrollcommand=self.scrollbar_area_vertical.set,
                                 xscrollcommand=self.scrollbar_area_horizontal.set,
                                 selectmode='multiple',
                                 exportselection=0)
        self.scrollbar_area_vertical.config(command=self.list_area.yview)
        self.scrollbar_area_vertical.pack(side=RIGHT, fill=Y)
        self.scrollbar_area_horizontal.config(command=self.list_area.xview)
        self.scrollbar_area_horizontal.pack(side=BOTTOM, fill=X)
        self.list_area.pack(side=LEFT, fill=BOTH, expand=1)

        #### (3.2.3) Interface Residues Listbox
        self.residues_frame = Frame(descriptors.interior())
        self.residues_frame.pack(side=LEFT, padx=10)
        self.scrollbar_residues_vertical = Scrollbar(self.residues_frame,
                                                     orient="vertical")
        self.scrollbar_residues_horizontal = Scrollbar(self.residues_frame,
                                                       orient="horizontal")
        Label(self.residues_frame,
              text="Interface Residues").pack()
        self.list_residues = Listbox(self.residues_frame,
                                     height=18,
                                     width=25,
                                     yscrollcommand=self.scrollbar_residues_vertical.set,
                                     xscrollcommand=self.scrollbar_residues_horizontal.set,
                                     selectmode='multiple',
                                     exportselection=0)
        self.scrollbar_residues_vertical.config(command=self.list_residues.yview)
        self.scrollbar_residues_vertical.pack(side=RIGHT, fill=Y)
        self.scrollbar_residues_vertical.config(command=self.list_residues.xview)
        self.scrollbar_residues_horizontal.pack(side=BOTTOM, fill=X)
        self.list_residues.pack(side=LEFT, fill=BOTH, expand=1)



        #### (3.2.4) Listboxes Function Binding
        self.list_volume.bind('<ButtonRelease-1>',
                              lambda event,
                                     arg=self.list_volume,
                                     arg2=self.list_area:
                              self.show_cavities(event, arg, arg2))
        self.list_area.bind('<ButtonRelease-1>',
                            lambda event,
                                   arg=self.list_area,
                                   arg2=self.list_volume:
                            self.show_cavities(event, arg, arg2))
        self.list_residues.bind('<ButtonRelease-1>',
                        lambda event,
                               arg=self.list_residues:
                        self.show_residues(event, arg))

        ################################################################################################################
        #    (4) About Tab                                                                                             #
        ################################################################################################################

        # Set up About tab
        about = self.notebook.add('About')

        #### (4.1) About Frame #########################################################################################
        about_frame = Frame(about)
        about_frame.pack(expand=1, fill=BOTH)

        #### (4.1.1) Text Holder
        scrollbar_about = Scrollbar(about_frame,
                                    orient="vertical")
        text_holder = Text(about_frame,
                           yscrollcommand=scrollbar_about.set,
                           background="#ddddff")
        scrollbar_about.config(command=text_holder.yview)
        scrollbar_about.pack(side=RIGHT, fill=Y)

        # Define text
        text = u"""PyMOL parKVFinder Tools integrates PyMOL (http://PyMOL.org/) with parKVFinder <website>.

In the simplest case to run parKVFinder:

1) Load a target biomolecular structure into PyMOL.
2) Start PyMOL parKVFinder Tools plugin.
3) Ensure that parKVFinder executable path is correct on the "Program Locations" tab.
4) Click the "Run parKVFinder" button.

parKVFinder and PyMOL parKVFinder Tools was developed by:
- Jo{}o Victor da Silva Guerra
- Helder Veras Ribeiro Filho
- Leandro Oliveira Bortot
- Rodrigo Vargas Honorato
- Jos{} Geraldo de Carvalho Pereira
- Paulo Sergio Lopes de Oliveira (paulo.oliveira@lnbio.cnpem.br)

Brazilian Center for Research in Energy and Materials - CNPEM
Brazilian Biosciences National Laboratory - LNBio

Please refer and cite the original paper if you use it in a publication.

Citation:
    Jo{}o Victor da Silva Guerra, Helder Veras Ribeiro Filho, Leandro Oliveira Bortot, Rodrigo Vargas Honorato, Jos{} Geraldo de Carvalho Pereira, Paulo Sergio Lopes de Oliveira, ParKVFinder: A thread-level parallel approach in biomolecular cavity detection, SoftwareX, 2020, https://doi.org/10.1016/j.softx.2020.100606

Citation for PyMOL may be found here:
    http://pymol.sourceforge.net/faq.html#CITE
""".format(chr(0x00e3),
           chr(0x00e9),
           chr(0x00e3),
           chr(0x00e9))
        # Insert text
        text_holder.insert(END, text)
        text_holder.pack(side=LEFT, expand=1, fill=BOTH)
        text_holder.configure(state='disabled')

        # Show plugin window
        self.notebook.setnaturalsize()
        self.show_app_modal()

        # Refresh app
        self.refresh_listbox(self.pdb_file)

    ####################################################################################################################
    #                                                                                                                  #
    #    (5) Inner Class: Box Adjustment Frame Class                                                                   #
    #                                                                                                                  #
    ####################################################################################################################
    class BoxAdjustmentWidget(Frame):
        """
        Box Adjustment Frame Graphical Interface and Methods
        """
        #### (5.1) Box Adjustment Frame Attributes #####################################################################
        def __init__(self, parent):
            #### Create frame
            Frame.__init__(self, parent)
            self.parent = parent
            self.widget()

            #### Set default box attributes
            self.box_name = "box"
            self.x = 0.0
            self.y = 0.0
            self.z = 0.0
            self.lim1 = 0.0
            self.lim2 = 0.0
            self.lim3 = 0.0
            self.lim4 = 0.0
            self.lim5 = 0.0
            self.lim6 = 0.0
            self.angle1 = 0.0
            self.angle2 = 0.0
            self.padding = 3.5

        #### (5.2) Buttons Functions ###################################################################################
        def set_box(self, padding):
            """
            Create box coordinates, enable 'Delete Box' and 'Redraw Box' buttons and call draw_box function.
            :param padding: box padding value.
            """
            # Delete Box object in PyMOL
            if "box" in cmd.get_names("selections"):
                cmd.delete(self.box_name)
            # Get dimensions of selected residues
            selection = "sele"
            if selection in cmd.get_names("selections"):
                ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)
            else:
                ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent("")

            # Get center of each dimension (x, y, z)
            self.x = (minX + maxX)/2
            self.y = (minY + maxY)/2
            self.z = (minZ + maxZ)/2

            # Set Box variables
            self.lim1 = round(self.x - (minX - float(padding.get())), 1)
            self.lim2 = round((maxX + float(padding.get())) - self.x, 1)
            self.lim3 = round(self.y - (minY - float(padding.get())), 1)
            self.lim4 = round((maxY + float(padding.get())) - self.y, 1)
            self.lim5 = round(self.z - (minZ - float(padding.get())), 1)
            self.lim6 = round((maxZ + float(padding.get())) - self.z, 1)
            self.angle1 = 0
            self.angle2 = 0
            self.padding = float(padding.get())

            # Setting values and enabling scrollbar
            self.min_x.setvalue(self.lim1)
            self.max_x.setvalue(self.lim2)
            self.min_y.setvalue(self.lim3)
            self.max_y.setvalue(self.lim4)
            self.min_z.setvalue(self.lim5)
            self.max_z.setvalue(self.lim6)
            self.ang_s.setvalue(self.angle1)
            self.ang2_s.setvalue(self.angle2)

            # Draw box
            self.draw_box()
            self.dbbutton.configure(state='disabled')
            self.rbbutton.configure(state='normal')

        def delete_box(self):
            """
            Delete box object, disable 'Delete Box' and 'Redraw Box' buttons and enable 'Draw Box' button.
            """
            # Reset all box variables
            self.x = 0
            self.y = 0
            self.z = 0
            self.angle1 = 0
            self.angle2 = 0
            self.lim1 = 0
            self.lim2 = 0
            self.lim3 = 0
            self.lim4 = 0
            self.lim5 = 0
            self.lim6 = 0

            # Delete Box and Vertices objects in PyMOL
            cmd.delete("vertices")
            cmd.delete(self.box_name)

            # Set Box variables in the interface
            self.min_x.setvalue(self.lim1)
            self.max_x.setvalue(self.lim2)
            self.min_y.setvalue(self.lim3)
            self.max_y.setvalue(self.lim4)
            self.min_z.setvalue(self.lim5)
            self.max_z.setvalue(self.lim6)
            self.ang_s.setvalue(self.angle1)
            self.ang2_s.setvalue(self.angle2)

            # Change state of buttons in the interface
            self.dbbutton.configure(state='normal')
            self.rbbutton.configure(state='disabled')

        def draw_box(self):
            """
            Draw box in PyMOL interface.
            :return: box object.
            """

            # Convert angle
            self.angle1 = (float(self.angle1)/180.0)*3.1415926
            self.angle2 = (float(self.angle2)/180.0)*3.1415926

            # Get positions of box vertices
            z1 = self.lim1*math.sin(float(self.angle2)) \
                + self.lim3*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                - self.lim5*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y1 = -self.lim3*math.cos(float(self.angle1)) \
                + (-self.lim5)*math.sin(float(self.angle1)) \
                + self.y
            x1 = -self.lim1*math.cos(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z2 = (-self.lim2)*math.sin(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y2 = (-self.lim3)*math.cos(float(self.angle1)) \
                + (-self.lim5)*math.sin(float(self.angle1)) \
                + self.y
            x2 = self.lim2*math.cos(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z3 = -(-self.lim1)*math.sin(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y3 = self.lim4*math.cos(float(self.angle1)) \
                + (-self.lim5)*math.sin(float(self.angle1)) \
                + self.y
            x3 = (-self.lim1)*math.cos(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z4 = -(-self.lim1)*math.sin(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y4 = (-self.lim3)*math.cos(float(self.angle1)) \
                + self.lim6*math.sin(float(self.angle1)) \
                + self.y
            x4 = (-self.lim1)*math.cos(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z5 = (-self.lim2)*math.sin(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y5 = self.lim4*math.cos(float(self.angle1)) \
                + (-self.lim5)*math.sin(float(self.angle1)) \
                + self.y
            x5 = self.lim2*math.cos(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + (-self.lim5)*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z6 = (-self.lim2)*math.sin(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y6 = (-self.lim3)*math.cos(float(self.angle1)) \
                + self.lim6*math.sin(float(self.angle1)) \
                + self.y
            x6 = self.lim2*math.cos(float(self.angle2)) \
                - (-self.lim3)*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z7 = -(-self.lim1)*math.sin(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y7 = self.lim4*math.cos(float(self.angle1)) \
                + self.lim6*math.sin(float(self.angle1)) \
                + self.y
            x7 = (-self.lim1)*math.cos(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x
            z8 = (-self.lim2)*math.sin(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.cos(float(self.angle2)) \
                + self.z
            y8 = self.lim4*math.cos(float(self.angle1)) \
                + self.lim6*math.sin(float(self.angle1)) \
                + self.y
            x8 = self.lim2*math.cos(float(self.angle2)) \
                - self.lim4*math.sin(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.lim6*math.cos(float(self.angle1))*math.sin(float(self.angle2)) \
                + self.x

            # Create box object
            pymol.stored.list = []
            if self.box_name in cmd.get_names("selections"):
                cmd.iterate(self.box_name, "stored.list.append((name, color))", quiet=1)
            list_color = pymol.stored.list
            cmd.delete(self.box_name)
            if len(list_color) > 0:
                for item in list_color:
                    at_name = item[0]
                    at_c = item[1]
                    cmd.set_color(at_name+"color", cmd.get_color_tuple(at_c))
            else:
                for at_name in ["v2", "v3", "v4", "v5", "v6", "v7", "v8", "v1x", "v1y", "v1z", "v2x", "v3y", "v4z"]:
                    cmd.set_color(at_name+"color", [0.86, 0.86, 0.86])

            # Create vertices
            cmd.pseudoatom(self.box_name, name="v2", pos=[x2, y2, z2], color="v2color")
            cmd.pseudoatom(self.box_name, name="v3", pos=[x3, y3, z3], color="v3color")
            cmd.pseudoatom(self.box_name, name="v4", pos=[x4, y4, z4], color="v4color")
            cmd.pseudoatom(self.box_name, name="v5", pos=[x5, y5, z5], color="v5color")
            cmd.pseudoatom(self.box_name, name="v6", pos=[x6, y6, z6], color="v6color")
            cmd.pseudoatom(self.box_name, name="v7", pos=[x7, y7, z7], color="v7color")
            cmd.pseudoatom(self.box_name, name="v8", pos=[x8, y8, z8], color="v8color")

            # Connect vertices
            cmd.select("vertices", "(name v3,v7)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v2,v6)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v5,v8)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v2,v5)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v4,v6)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v4,v7)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v3,v5)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v6,v8)")
            cmd.bond("vertices", "vertices")
            cmd.select("vertices", "(name v7,v8)")
            cmd.bond("vertices", "vertices")
            cmd.pseudoatom(self.box_name, name="v1x", pos=[x1, y1, z1], color='red')
            cmd.pseudoatom(self.box_name, name="v2x", pos=[x2, y2, z2], color='red')
            cmd.select("vertices", "(name v1x,v2x)")
            cmd.bond("vertices", "vertices")
            cmd.pseudoatom(self.box_name, name="v1y", pos=[x1, y1, z1], color='forest')
            cmd.pseudoatom(self.box_name, name="v3y", pos=[x3, y3, z3], color='forest')
            cmd.select("vertices", "(name v1y,v3y)")
            cmd.bond("vertices", "vertices")
            cmd.pseudoatom(self.box_name, name="v4z", pos=[x4, y4, z4], color='blue')
            cmd.pseudoatom(self.box_name, name="v1z", pos=[x1, y1, z1], color='blue')
            cmd.select("vertices", "(name v1z,v4z)")
            cmd.bond("vertices", "vertices")
            cmd.delete("vertices")

        def redraw_box(self, padding):
            """
            Redraw box in PyMOL interface.
            :param padding: box padding.
            :return: box object.
            """

            if "sele" in cmd.get_names("selections"):
                # Get dimensions of selected residues
                selection = "sele"
                ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)

                if float(self.min_x.get()) != self.lim1 or float(self.max_x.get()) != self.lim2 or float(self.min_y.get()) != self.lim3 or float(self.max_y.get()) != self.lim4 or float(self.min_z.get()) != self.lim5 or float(self.max_z.get()) != self.lim6 or float(self.ang_s.get()) != self.angle1 or float(self.ang2_s.get()) != self.angle2:
                    self.lim1 = float(self.min_x.get())
                    self.lim2 = float(self.max_x.get())
                    self.lim3 = float(self.min_y.get())
                    self.lim4 = float(self.max_y.get())
                    self.lim5 = float(self.min_z.get())
                    self.lim6 = float(self.max_z.get())
                    self.angle1 = float(self.ang_s.get())
                    self.angle2 = float(self.ang2_s.get())

                # Padding was altered and/or selection was altered
                else:
                    # Get center of each dimension (x, y, z)
                    self.x = (minX + maxX)/2
                    self.y = (minY + maxY)/2
                    self.z = (minZ + maxZ)/2

                    # Set Box variables
                    self.lim1 = round(self.x - (minX - float(padding.get())), 1) + (float(self.min_x.get()) - self.lim1)
                    self.lim2 = round((maxX + float(padding.get())) - self.x, 1) + (float(self.max_x.get()) - self.lim2)
                    self.lim3 = round(self.y - (minY - float(padding.get())), 1) + (float(self.min_y.get()) - self.lim3)
                    self.lim4 = round((maxY + float(padding.get())) - self.y, 1) + (float(self.max_y.get()) - self.lim4)
                    self.lim5 = round(self.z - (minZ - float(padding.get())), 1) + (float(self.min_z.get()) - self.lim5)
                    self.lim6 = round((maxZ + float(padding.get())) - self.z, 1) + (float(self.max_z.get()) - self.lim6)
                    self.angle1 = 0 + float(self.ang_s.get())
                    self.angle2 = 0 + float(self.ang2_s.get())
                    self.padding = float(padding.get())

                    # Setting values and enabling scrollbar
                    self.min_x.setvalue(self.lim1)
                    self.max_x.setvalue(self.lim2)
                    self.min_y.setvalue(self.lim3)
                    self.max_y.setvalue(self.lim4)
                    self.min_z.setvalue(self.lim5)
                    self.max_z.setvalue(self.lim6)
                    self.ang_s.setvalue(self.angle1)
                    self.ang2_s.setvalue(self.angle2)

            else:
                # Not considering selected residues
                if float(self.min_x.get()) != self.lim1 or float(self.max_x.get()) != self.lim2 or float(self.min_y.get()) != self.lim3 or float(self.max_y.get()) != self.lim4 or float(self.min_z.get()) != self.lim5 or float(self.max_z.get()) != self.lim6 or float(self.ang_s.get()) != self.angle1 or float(self.ang2_s.get()) != self.angle2:
                    self.lim1 = float(self.min_x.get())
                    self.lim2 = float(self.max_x.get())
                    self.lim3 = float(self.min_y.get())
                    self.lim4 = float(self.max_y.get())
                    self.lim5 = float(self.min_z.get())
                    self.lim6 = float(self.max_z.get())
                    self.angle1 = float(self.ang_s.get())
                    self.angle2 = float(self.ang2_s.get())

                # Padding was altered
                if self.padding != float(padding.get()):
                    # Prepare dimensions without old padding
                    minX = self.padding - float(self.min_x.get())
                    maxX = float(self.max_x.get()) - self.padding
                    minY = self.padding - float(self.min_y.get())
                    maxY = float(self.max_y.get()) - self.padding
                    minZ = self.padding - float(self.min_z.get())
                    maxZ = float(self.max_z.get()) - self.padding

                    # Get center of each dimension (x, y, z)
                    self.x = (minX + maxX)/2
                    self.y = (minY + maxY)/2
                    self.z = (minZ + maxZ)/2

                    # Set box variables
                    self.lim1 = round(self.x - (minX - float(padding.get())), 1)
                    self.lim2 = round((maxX + float(padding.get())) - self.x, 1)
                    self.lim3 = round(self.y - (minY - float(padding.get())), 1)
                    self.lim4 = round((maxY + float(padding.get())) - self.y, 1)
                    self.lim5 = round(self.z - (minZ - float(padding.get())), 1)
                    self.lim6 = round((maxZ + float(padding.get())) - self.z, 1)
                    self.angle1 = float(self.ang_s.get())
                    self.angle2 = float(self.ang2_s.get())
                    self.padding = float(padding.get())

                    # Set Box variables in the interface
                    self.min_x.setvalue(self.lim1)
                    self.max_x.setvalue(self.lim2)
                    self.min_y.setvalue(self.lim3)
                    self.max_y.setvalue(self.lim4)
                    self.min_z.setvalue(self.lim5)
                    self.max_z.setvalue(self.lim6)
                    self.ang_s.setvalue(self.angle1)
                    self.ang2_s.setvalue(self.angle2)

            # Redraw Box
            self.draw_box()

        def _custom_real_counter(self, text, factor, increment):
            """
            Custom counter function for Pmw.Counter increment and decrement buttons.
            :param text: Current numeric value in Pmw.EntryField.
            :param factor: Value of 1 or -1 to multiply the increment parameter.
            :param increment: Value to increment in text parameter.
            :return: Numeric value to show in Pmw.EntryField.
            """
            return "{0:.1f}".format(float(text)+increment*factor)

        def widget(self):
            """
            Create widget in Box Adjustment Frame of PyMOL parKVFinder Tools to interact with box object.
            """
            Label(self, text='Select residues and press Draw Box').pack()

            # Draw Box Button
            buttons = Frame(self)
            draw_box_button = lambda: self.set_box(padding)
            self.dbbutton = Button(buttons, text="Draw Box", command=draw_box_button)
            self.dbbutton.pack(side=LEFT)

            # Delete Box Button
            delete_boxbutton = lambda: self.delete_box()
            Button(buttons, text="Delete Box", command=delete_boxbutton).pack(side=LEFT)

            # Delete Box Button
            redraw_box_button = lambda: self.redraw_box(padding)
            self.rbbutton = Button(buttons, text="Redraw Box", command=redraw_box_button)
            self.rbbutton.pack(side=LEFT)
            self.rbbutton.configure(state='disabled')
            buttons.pack()

            # Padding
            padding = Pmw.Counter(self,
                                  labelpos='w',
                                  label_text=u"Padding ({}):".format(chr(0x212b)),
                                  datatype={'counter': self._custom_real_counter},
                                  entryfield_validate={'validator': 'real', 'min': 0, 'max': 10},
                                  increment=0.1,
                                  entry_width=5,
                                  entryfield_value=3.5)
            padding.pack(fill=X)

            # Minimum X
            self.min_x = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Minimum X ({}):".format(chr(0x212b)),
                                     label_foreground='red',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': -50, 'max': 50},
                                     increment=0.1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.min_x.pack(fill=X)

            # Maximum X
            self.max_x = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Maximum X ({}):".format(chr(0x212b)),
                                     label_foreground='red',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': -50, 'max': 50},
                                     increment=0.1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.max_x.pack(fill=X)

            # Minimum Y
            self.min_y = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Minimum Y ({}):".format(chr(0x212b)),
                                     label_foreground='green',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': -50, 'max': 50},
                                     increment=0.1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.min_y.pack(fill=X)

            # Maximum Y
            self.max_y = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Maximum Y ({}):".format(chr(0x212b)),
                                     label_foreground='green',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': -50, 'max': 50},
                                     increment=0.1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.max_y.pack(fill=X)

            # Minimum Z
            self.min_z = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Minimum Z ({}):".format(chr(0x212b)),
                                     label_foreground='blue',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': -50, 'max': 50},
                                     increment=0.1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.min_z.pack(fill=X)

            # Maximum Z
            self.max_z = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Maximum Z ({}):".format(chr(0x212b)),
                                     label_foreground='blue',
                                     datatype={'counter': self._custom_real_counter},
                                     entryfield_validate={'validator': 'real', 'min': -50, 'max': 50},
                                     increment=0.1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.max_z.pack(fill=X)

            # Angle 1
            self.ang_s = Pmw.Counter(self,
                                     labelpos='w',
                                     label_text=u"Angle 1 ({}):".format(chr(176)),
                                     label_foreground='black',
                                     datatype='real',
                                     entryfield_validate={'validator': 'real', 'min': -180, 'max': 180},
                                     increment=1,
                                     entry_width=5,
                                     entryfield_value=0)
            self.ang_s.pack(fill=X)

            # Angle 2
            self.ang2_s = Pmw.Counter(self,
                                      labelpos='w',
                                      label_text=u"Angle 2 ({}):".format(chr(176)),
                                      label_foreground='black',
                                      datatype='real',
                                      entryfield_validate={'validator': 'real', 'min': -180, 'max': 180},
                                      increment=1,
                                      entry_width=5,
                                      entryfield_value=0)
            self.ang2_s.pack(fill=X)

    ####################################################################################################################
    #                                                                                                                  #
    #    (6) Methods                                                                                                   #
    #                                                                                                                  #
    ####################################################################################################################

    #### (6.1) Custom counter function for Pmw.Counter #################################################################
    def _custom_real_counter(self, text, factor, increment):
        """
        Custom counter function for Pmw.Counter increment and decrement buttons.
        :param text: Current numeric value in Pmw.EntryField.
        :param factor: Value of 1 or -1 to multiply the increment parameter.
        :param increment: Value to increment in text parameter.
        :return: Numeric value to show in Pmw.EntryField.
        """
        return "{0:.1f}".format(float(text) + increment * factor)

    #### (6.2) Check Functions #########################################################################################
    def check_grid_method(self):
        """
        Check logical conditions of Grid Method values.
        :return: Set Resolution, Step Size and Grid Method to acceptable values.
        """
        if self.grid_method.get():
            self.resolution.setvalue('Off')
        else:
            self.step_size.setentry(0.0)
            if self.resolution.getcurselection() == 'Off':
                self.resolution.setvalue('Low')

    def check_ligand_mode(self):
        """
        Check logical conditions to enable or disable features of Ligand Adjustment Group.
        :return: Enable/disable buttons and erase load structures in PyMOL of Listbox.
        """
        if self.ligand_mode.get():
            self.ligand_pdb_refresh_button.configure(state='normal')
            self.ligand_pdb_upload_button.configure(state='normal')
        else:
            self.ligand_pdb_refresh_button.configure(state='disabled')
            self.ligand_pdb_upload_button.configure(state='disabled')
            self.ligand_file.delete(0, self.ligand_file.size())

    def check_step_size_status(self, value):
        """
        Check logical conditions to resolutions flags.
        :param value: Resolution flag set.
        :return: Set Resolution, Step Size and Grid Method to acceptable values.
        """
        if value != 'Off':
            self.step_size.setentry(0.0)
            self.grid_method.set(0)
        else:
            self.grid_method.set(1)

    #### (6.3) Buttons Functions #######################################################################################
    def upload_pdb(self, listbox):
        """
        Upload PDB file with tkFileDialog.
        :param listbox: Listbox to be refreshed.
        :return: Load a structure in PyMOL and show it in the listbox.
        """
        filename = tkFileDialog.askopenfilename(
            filetypes=[
                ("Ligand PDB File", "*.pdb"),
                ("Ligand CIF File", "*.cif"),
                ("All files", "*")
            ]
        )
        if len(filename) != 0:
            cmd.load(filename)
            self.refresh_listbox(listbox)

    def refresh_listbox(self, listbox):
        """
        Refresh Listbox with loaded structures in PyMOL.
        :param listbox: Listbox to be refreshed.
        :return: Show loaded structures in the listbox.
        """
        listbox.delete(0, listbox.size())
        for item in cmd.get_names("all"):
            if cmd.get_type(item) == "object:molecule" and \
                    item != "box" and \
                    item != "grid" and \
                    item != "cavities" and \
                    item != "residues" and \
                    item[-16:] != ".KVFinder.output" and \
                    item != 'target_exclusive':
                listbox.insert(listbox.size(), item)
        listbox.selection_set(0, 0)

    def select_file(self, file, title, file_types=[('All files', '*')]):
        """
        Open a tkFileDialog to select a file.
        :param file: Entry to show path of selected file.
        :param title: Title to show in the tkFileDialog header.
        :param file_types: Types of files to be show in tkFileDialog.
        :return: Fill Entry with selected file path.
        """
        file_path = tkFileDialog.askopenfilename(title=title,
                                                 parent=self.w,
                                                 filetypes=file_types)
        if len(file_path) != 0:
            file.configure(state='normal')
            file.delete(0, len(file.get()))
            file.insert(0, file_path)
            file.configure(state='readonly')

    def select_directory(self, directory, title):
        """
        Open a tkFileDialog to select a directory.
        :param directory: Entry to show path of selected directory.
        :param title: Title to show in the tkFileDialog header.
        :return: Fill Entry with selected directory path.
        """
        directory_path = tkFileDialog.askdirectory(mustexist=1,
                                                   title=title,
                                                   parent=self.w)
        if len(directory_path) != 0:
            directory.configure(state='normal')
            directory.delete(0, len(directory.get()))
            directory.insert(0, directory_path)
            directory.configure(state='readonly')

    #### (6.4) Handle Box Adjustment Frame Function ####################################################################
    def update_box_adjustment_frame(self, event):
        """
        Handle the Box Adjustment Frame based on the selected Search Procedure Method.
        :param event: Click in Search Procedure Radiobutton.
        :return: Create or destroy Box Adjustment Frame.
        """
        if event == 'Whole Protein':
            cmd.delete("box")
            if self.box_adjustment_frame.winfo_exists():
                self.box_adjustment_frame.destroy()
                self.box_adjustment_frame = Frame(self.adjustment_method_group_frame1, width=350)
                self.box_adjustment_frame.pack(padx=10, pady=5, fill=X)
        else:
            if self.box_adjustment_frame.winfo_exists() and \
                    len(self.box_adjustment_frame.winfo_children()) == 0:
                # Create Box Adjustment Group

                self.ba_group = Pmw.Group(self.box_adjustment_frame,
                                          tag_text='Box Adjustment')
                self.ba_group.pack(padx=15, pady=5, fill=X)
                # Create Box adjustment
                self.box_adjustment_widget = self.BoxAdjustmentWidget(self.ba_group.interior())
                self.box_adjustment_widget.pack(padx=10, pady=5, fill=X)

    #### (6.5) Results Visualization Functions #########################################################################
    def clean_results(self):
        """
        Erase Information Group and Descriptor Group in Results Visualization Tab.
        :return: None
        """
        #### Clean Information Group
        ##  Erase Input File
        self.inp_file.configure(
            text="Input File: "
        )
        ## Erase Ligand File
        self.lig_file.configure(
            text="Ligand File: "
        )
        ## Erase Output File
        self.output_file.configure(
            text="Output File: "
        )
        ## Erase Step Size
        self.step_size_info.configure(
            text=u"Step Size ({}): ".format(
                chr(0x212b)
            )
        )

        #### Clean Descriptor Group
        ## Erase Volume Listbox
        self.list_volume.delete(
            0,
            self.list_volume.size()
        )
        ## Erase Area Listbox
        self.list_area.delete(
            0,
            self.list_area.size()
        )
        ## Erase Interface Residues Listbox
        self.list_residues.delete(
            0,
            self.list_residues.size()
        )

    def refresh_residues_information(self):
        """
        Refresh Residues Listbox with residues that form each cavity contained in <output>.KVFinder.results.toml.
        :return: None
        """
        # Insert cavity index with residues information in Residues Listbox
        for cavity_index in sorted(results['RESULTS']['RESIDUES'].keys()):
            self.list_residues.insert(self.list_residues.size(), cavity_index.rstrip('\n'))

    def refresh_spatial_information(self):
        """
        Refresh Volume and Area Listbox for each cavity contained in <output>.KVFinder.results.toml.
        :return: None
        """
        #### Get all cavity indexes
        cavity_indexes = sorted(results['RESULTS']['VOLUME'].keys())
        for cavity_index in cavity_indexes:
            # Include Volume
            self.list_volume.insert(
                self.list_volume.size(),
                "{}: {}".format(cavity_index, results['RESULTS']['VOLUME'][cavity_index])
            )
            # Include Area
            self.list_area.insert(
                self.list_area.size(),
                "{}: {}".format(cavity_index, results['RESULTS']['AREA'][cavity_index])
            )

    def load_results(self, loaded_results):
        """
        Load cavity detection and characterization results from <output>.KVFinder.results.toml file.
        :param self.results_file: Path to <output>.KVFinder.results.toml file to be loaded.
        :type self.results_file: Tkinter.Entry object
        :returns: None
        """
        #### Check <output>.KVFinder.results.toml file exists and has appropriate extension
        if os.path.exists(self.results_file.get()) and self.results_file.get().endswith('KVFinder.results.toml'):
            sys.stdout.write("Loading results from: {}\n\n".format(self.results_file.get()))
        else:
            tkMessageBox.showerror("Error",
                                   "Results file cannot be opened! Check results file path.",
                                   parent=self.w
                                   )
            return False

        #### Create a global variable for results
        global results

        #### Read <output>.KVFinder.results.toml file
        results = toml.load(self.results_file.get())

        #### Clean results
        self.clean_results()

        #### Update Information Group
        ## Input PDB File
        self.inp_file.configure(
            text="Input File: {}".format(
                results['FILES_PATH']['INPUT']
            )
        )
        ## Ligand PDB File
        self.lig_file.configure(
            text="Ligand File: {}".format(
                results['FILES_PATH']['LIGAND']
            )
        )
        ## Output PDB File
        self.output_file.configure(
            text="Output file: {}".format(
                results['FILES_PATH']['OUTPUT']
            )
        )
        ## Step Size
        self.step_size_info.configure(
            text=u"Step Size ({}): {:.2f}".format(
                chr(0x212b),
                results['PARAMETERS']['STEP_SIZE']
            )
        )

        #### Update objects
        ## Delete old results
        cmd.delete("cavities")
        cmd.delete("residues")
        cmd.frame(1)
        ## Remove previous results in objects with same name
        for x in cmd.get_names("all"):
            if results['FILES_PATH']['OUTPUT'].split('/')[-1].replace('.pdb', '') == x:
                cmd.delete(x)
        ## Load <output>.KVFinder.output.pdb
        if os.path.exists(results['FILES_PATH']['OUTPUT']):
            # cmd.set("auto_zoom", 0)
            cmd.load(results['FILES_PATH']['OUTPUT'],
                     results['FILES_PATH']['OUTPUT'].split('/')[-1].replace('.pdb', ''),
                     zoom=0)
            cmd.hide('everything', results['FILES_PATH']['OUTPUT'].split('/')[-1].replace('.pdb', ''))
            cmd.show('nonbonded', results['FILES_PATH']['OUTPUT'].split('/')[-1].replace('.pdb', ''))
            # cmd.set("auto_zoom", 1)

        if loaded_results:
            cmd.disable()
            cmd.enable(results['FILES_PATH']['OUTPUT'].split('/')[-1].replace('.pdb', ''))
            self.last_input.set(results['FILES_PATH']['INPUT'].split('/')[-1].replace('.pdb', ''),)
            for x in cmd.get_names("all"):
                if results['FILES_PATH']['INPUT'].split('/')[-1].replace('.pdb', '') == x:
                    cmd.delete(x)
            cmd.load(results['FILES_PATH']['INPUT'],
                     results['FILES_PATH']['INPUT'].split('/')[-1].replace('.pdb', ''),
                     zoom=0)
            self.last_object.set(results['FILES_PATH']['OUTPUT'].split('/')[-1].replace('.pdb', ''))
            if results['FILES_PATH']['LIGAND'] != '-':
                for x in cmd.get_names("all"):
                    if results['FILES_PATH']['LIGAND'].split('/')[-1].replace('.pdb', '') == x:
                        cmd.delete(x)
                cmd.load(results['FILES_PATH']['LIGAND'],
                         results['FILES_PATH']['LIGAND'].split('/')[-1].replace('.pdb', ''),
                         zoom=0)
                self.last_ligand.set(results['FILES_PATH']['LIGAND'].split('/')[-1].replace('.pdb', ''))
            self.loaded_results.set(1)
        else:
            self.loaded_results.set(0)

        # Load spatial information
        self.refresh_spatial_information()

        # Load residues information
        self.refresh_residues_information()

    def show_residues(self, event, arg):
        """
        Show residues of selected cavities in Interface Residues Listbox.
        :param event: Click in Interface Residues Listbox.
        :param arg: Listbox of Interface Residues.
        :return: residues object in PyMOL with residues of selected cavities.
        """
        residues_list = []
        for item in arg.curselection():
            for i in results['RESULTS']['RESIDUES'][arg.get(item)]:
                if i not in residues_list:
                    residues_list.append(i)

        cmd.set("auto_zoom", 0)
        cmd.delete("res")
        cmd.delete("residues")
        if len(arg.curselection()) < 1:
            return

        control = 0
        for item in cmd.get_names("all"):
            if item == self.last_object.get():
                control = 1
        if control == 0:
            return

        # Select residues around chosen cavities
        # Create command
        command = StringVar()
        command.set("{} and".format(self.last_input.get()))
        while len(residues_list) > 1:
            command.set("{} (resid {} and chain {}) or".format(command.get(),
                                                               residues_list[0][0],
                                                               residues_list[0][1]))
            residues_list.remove(residues_list[0])
        command.set("{} (resid {} and chain {})".format(command.get(),
                                                        residues_list[0][0],
                                                        residues_list[0][1]))
        residues_list.remove(residues_list[0])

        # Create residues object with residues around chosen cavities
        cmd.select("res", command.get())
        cmd.create("residues", "res")
        cmd.delete("res")

        # Highlight residues
        cmd.hide("everything", "residues")
        cmd.show('sticks', 'residues')
        cmd.disable(self.last_object.get())
        cmd.enable(self.last_object.get())
        cmd.set("auto_zoom", 1)

    def show_cavities(self, event, arg, arg2):
        """
        Show cavity points of selected cavities in Volume or Area Listbox.
        :param event: Click in Volume or Area Listbox.
        :param arg: Listbox of Volume or Area.
        :param arg2: Listbox of Volume or Area.
        :return: cavities object in PyMOL with cavity points of selected cavities.
        """
        cavs = []
        for item in arg.curselection():
            cavs.append(arg.get(item)[0:3])
            arg2.selection_set(item, item)

        for item in arg2.curselection():
            if item not in arg.curselection():
                arg2.selection_clear(item, item)

        cmd.set("auto_zoom", 0)
        cmd.delete("cavs")  # cmd.delete("cavities")
        cmd.delete("cavities")  # cmd.delete("cavities.KVFinder.output")
        if len(arg.curselection()) < 1:
            return

        control = 0
        for item in cmd.get_names("all"):
            if item == self.last_object.get():
                control = 1
        if control == 0:
            return

        # Color filling cavity points as blue nonbounded
        command = StringVar()
        while len(cavs) > 0:
            command.set("{} and (resname ".format(self.last_object.get()))
            for item in cavs[:50]:
                command.set(command.get()+item+',')
                cavs.remove(item)
            command.set(command.get()[:-1]+')')
            cmd.select("cavs", command.get())

        # Create cavities object as blue nonbounded
        cmd.create("cavities", "cavs")  # cmd.create("cavities.KVFinder.output", "cavities")
        cmd.delete("cavs")  # cmd.delete("cavities")
        cmd.color("blue", "cavities")  # cmd.color("blue", "cavities.KVFinder.output")
        cmd.show('nonbonded', 'cavities')  # cmd.show('nonbonded', 'cavities.KVFinder.output')

        # Color surface cavity points as red nb_spheres
        cmd.select("cavs", "cavities and name HS")
        cmd.color("red", "cavs")
        cmd.show("nb_spheres", "cavs")
        cmd.delete("cavs")

        # Reset cavities output object
        cmd.disable(self.last_object.get())
        cmd.enable(self.last_object.get())
        cmd.set("auto_zoom", 1)

    #### (6.6) Dialog Buttons Functions ################################################################################
    def create_parKVFinder_input(self):
        """
        Check logical conditions to create parameters.toml of set parameters.
        :return: A parameters.toml file with set parameters in PyMOL parKVFinder Tools.
        """
        # Conditions tests
        # PDB file
        if not self.pdb_file.curselection():
            tkMessageBox.showerror("Error",
                                   "Select a input file",
                                   parent=self.parent)
            return
        # Ligand file
        if self.ligand_mode.get():
            if not self.ligand_file.curselection():
                tkMessageBox.showerror("Error",
                                       "Select a ligand file",
                                       parent=self.parent)
                return

        # Create output folder
        # Check if KV_Files directory exists. If not, create KV_Files directory in provided path
        if os.path.isdir('{}/KV_Files'.format(self.output.get())) == 0:
            subprocess.call(['mkdir', '{}/KV_Files'.format(self.output.get())])
        if os.path.isdir('{}/KV_Files/{}'.format(self.output.get(), self.base_name.getvalue())) == 0:
            subprocess.call(['mkdir', '{}/KV_Files/{}'.format(self.output.get(), self.base_name.getvalue())])

        # Create pdb file
        for x in cmd.get_names("all"):
            if x == self.pdb_file.get(self.pdb_file.curselection()):
                pdb_name = "{}/KV_Files/{}/{}.KVFinder.input.pdb".format(self.output.get(),
                                                                         self.base_name.getvalue(),
                                                                         self.base_name.getvalue())
                cmd.save(pdb_name, self.pdb_file.get(self.pdb_file.curselection()), 0, "pdb")

        # Create ligand file
        if self.ligand_mode.get():
            for x in cmd.get_names("all"):
                if x == self.pdb_file.get(self.pdb_file.curselection()):
                    ligand_name = "{}/KV_Files/{}/ligand_{}.KVFinder.input.pdb".format(self.output.get(),
                                                                                    self.base_name.getvalue(),
                                                                                    self.base_name.getvalue())
                    cmd.save(ligand_name, self.ligand_file.get(self.ligand_file.curselection()), 0, "pdb")
        else:
            ligand_name = "-"

        # Create parameters
        param = u"""# TOML configuration file for parKVFinder software.

title = "parKVFinder parameters file"

[FILES_PATH]
# The path of van der Waals radii dictionary for parKVFinder.
dictionary = \"{}\"
# The path of the input PDB file.
pdb = \"{}\"
# The path of the output directory.
output = \"{}\"
# Base name for output files.
base_name = \"{}\"
# The path for the ligand's PDB file.
ligand = \"{}\"

[SETTINGS]
# Settings for parKVFinder software

\t[SETTINGS.modes]
\t# Whole protein mode defines the search space as the whole protein.
\twhole_protein_mode = {}
\t# Box adjustment mode defines the search space as the box drawn in PyMOL.
\tbox_mode = {}
\t# Resolution mode implicitly sets the step size (grid spacing) of the 3D grid.
\t# If set to High, sets a voxel volume of 0.2. If set to Medium, sets a voxel volume of 0.1. \
If set to Low, it sets a voxel volume of 0.01. If set to Off, the step size must be set explicitly.
\tresolution_mode = \"{}\"
\t# Surface mode defines the type of surface representation to be applied, \
van der Waals molecular surface (true) or solvent accessible surface (false).
\tsurface_mode = {}
\t# Cavity output mode defines whether cavities are exported to the output PDB file as filled cavities \
(true) or filtered cavities (false).
\tkvp_mode = {}
\t# Ligand adjustment mode defines the search space around the ligand.
\tligand_mode = {}

\t[SETTINGS.step_size]
\t# Sets the 3D grid spacing. It directly affects accuracy and runtime.
\tstep_size = {}

\t[SETTINGS.probes]
\t# parKVFinder works with a dual probe system. \
A smaller probe, called Probe In, and a bigger one, called Probe Out, rolls around the protein.
\t# Points reached by the Probe In, but not the Probe Out are considered cavity points.
\t# Sets the Probe In diameter. Default: 1.4 angstroms.
\tprobe_in = {}
\t# Sets the Probe Out diameter. Default: 4.0 angstroms.
\tprobe_out = {}

\t[SETTINGS.cutoffs]
\t# Sets a volume cutoff for the detected cavities. Default: 5.0 angstroms.
\tvolume_cutoff = {}
\t# Sets a distance cutoff for a search space around the ligand in ligand adjustment mode. Default: 5.0 angstroms.
\tligand_cutoff = {}
\t# Sets a removal distance for the cavity frontier, which is defined by comparing Probe In and Probe Out surfaces. \
Default: 2.4 angstroms.
\tremoval_distance = {}
""".format(self.dictionary.get(),
           pdb_name,
           self.output.get(),
           self.base_name.getvalue(),
           ligand_name,
           "true" if self.search_procedure_options[self.search_procedure.getvalue()] else "false",
           "false" if self.search_procedure_options[self.search_procedure.getvalue()] else "true",
           self.resolution.getcurselection(),
           "true" if self.surface_representation[self.surface.getvalue()] else "false",
           "true" if self.cavity_representation_options[self.cavity_representation.getvalue()] else "false",
           "true" if self.ligand_mode.get() else "false",
           self.step_size.get(),
           self.probe_in.get(),
           self.probe_out.get(),
           self.volume_cutoff.get(),
           self.ligand_cutoff.get(),
           self.removal_distance.get())

        # Box adjustment mode is activated, so get box parameters
        if len(self.box_adjustment_frame.winfo_children()):
            # print("Box_mode ON")
            lim1 = float(self.box_adjustment_widget.lim1)
            lim2 = float(self.box_adjustment_widget.lim2)
            lim3 = float(self.box_adjustment_widget.lim3)
            lim4 = float(self.box_adjustment_widget.lim4)
            lim5 = float(self.box_adjustment_widget.lim5)
            lim6 = float(self.box_adjustment_widget.lim6)
            angle1 = float(self.box_adjustment_widget.angle1)
            angle2 = float(self.box_adjustment_widget.angle2)
            x = float(self.box_adjustment_widget.x)
            y = float(self.box_adjustment_widget.y)
            z = float(self.box_adjustment_widget.z)
        # Whole protein mode is activated, so create default box parameters
        else:
            # print("Whole protein ON")
            lim1 = 0.0
            lim2 = 0.0
            lim3 = 0.0
            lim4 = 0.0
            lim5 = 0.0
            lim6 = 0.0
            angle1 = 0.0
            angle2 = 0.0
            x = 0.0
            y = 0.0
            z = 0.0

        # Prepare points (bP1, bP2, bP3, bP4, P1, P2, P3, P4)
        # bP1
        bZ1 = lim1*math.sin(angle2) \
            + lim3*math.sin(angle1)*math.cos(angle2) \
            - lim5*math.cos(angle1)*math.cos(angle2) \
            + z
        bY1 = (-lim3)*math.cos(angle1) \
            + (-lim5)*math.sin(angle1) \
            + y
        bX1 = (-lim1)*math.cos(angle2) \
            - (-lim3)*math.sin(angle1)*math.sin(angle2) \
            + (-lim5)*math.cos(angle1)*math.sin(angle2) \
            + x
        # bP2
        bZ2 = (-lim2)*math.sin(angle2) \
            - (-lim3)*math.sin(angle1)*math.cos(angle2) \
            + (-lim5)*math.cos(angle1)*math.cos(angle2) \
            + z
        bY2 = (-lim3)*math.cos(angle1) \
            + (-lim5)*math.sin(angle1) \
            + y
        bX2 = lim2*math.cos(angle2) \
            - (-lim3)*math.sin(angle1)*math.sin(angle2) \
            + (-lim5)*math.cos(angle1)*math.sin(angle2) \
            + x
        # bP3
        bZ3 = -(-lim1)*math.sin(angle2) \
            - lim4*math.sin(angle1)*math.cos(angle2) \
            + (-lim5)*math.cos(angle1)*math.cos(angle2) \
            + z
        bY3 = lim4*math.cos(angle1) \
            + (-lim5)*math.sin(angle1) \
            + y
        bX3 = (-lim1)*math.cos(angle2) \
            - lim4*math.sin(angle1)*math.sin(angle2) \
            + (-lim5)*math.cos(angle1)*math.sin(angle2) \
            + x
        # bP4
        bZ4 = -(-lim1)*math.sin(angle2) \
            - (-lim3)*math.sin(angle1)*math.cos(angle2) \
            + lim6*math.cos(angle1)*math.cos(angle2) \
            + z
        bY4 = (-lim3)*math.cos(angle1) \
            + lim6*math.sin(angle1) \
            + y
        bX4 = (-lim1)*math.cos(angle2) \
            - (-lim3)*math.sin(angle1)*math.sin(angle2) \
            + lim6*math.cos(angle1)*math.sin(angle2) \
            + x
        # P1
        Z1 = (lim1+float(self.probe_out.get()))*math.sin(angle2) \
            + (lim3+float(self.probe_out.get()))*math.sin(angle1)*math.cos(angle2) \
            - (lim5+float(self.probe_out.get()))*math.cos(angle1)*math.cos(angle2) \
            + z
        Y1 = (-lim3-float(self.probe_out.get()))*math.cos(angle1) \
            + (-lim5-float(self.probe_out.get()))*math.sin(angle1) \
            + y
        X1 = (-lim1-float(self.probe_out.get()))*math.cos(angle2) \
            - (-lim3-float(self.probe_out.get()))*math.sin(angle1)*math.sin(angle2) \
            + (-lim5-float(self.probe_out.get()))*math.cos(angle1)*math.sin(angle2) \
            + x
        # P2
        Z2 = -(lim2+float(self.probe_out.get()))*math.sin(angle2) \
            - (-lim3-float(self.probe_out.get()))*math.sin(angle1)*math.cos(angle2) \
            + (-lim5-float(self.probe_out.get()))*math.cos(angle1)*math.cos(angle2) \
            + z
        Y2 = (-lim3-float(self.probe_out.get()))*math.cos(angle1) \
            + (-lim5-float(self.probe_out.get()))*math.sin(angle1) \
            + y
        X2 = (lim2+float(self.probe_out.get()))*math.cos(angle2) \
            - (-lim3-float(self.probe_out.get()))*math.sin(angle1)*math.sin(angle2) \
            + (-lim5-float(self.probe_out.get()))*math.cos(angle1)*math.sin(angle2) \
            + x
        # P3
        Z3 = -(-lim1-float(self.probe_out.get()))*math.sin(angle2) \
            - (lim4+float(self.probe_out.get()))*math.sin(angle1)*math.cos(angle2) \
            + (-lim5-float(self.probe_out.get()))*math.cos(angle1)*math.cos(angle2) \
            + z
        Y3 = (lim4+float(self.probe_out.get()))*math.cos(angle1) \
            + (-lim5-float(self.probe_out.get()))*math.sin(angle1) \
            + y
        X3 = (-lim1-float(self.probe_out.get()))*math.cos(angle2) \
            - (lim4+float(self.probe_out.get()))*math.sin(angle1)*math.sin(angle2) \
            + (-lim5-float(self.probe_out.get()))*math.cos(angle1)*math.sin(angle2) \
            + x
        # P4
        Z4 = -(-lim1-float(self.probe_out.get()))*math.sin(angle2) \
            - (-lim3-float(self.probe_out.get()))*math.sin(angle1)*math.cos(angle2) \
            + (lim6+float(self.probe_out.get()))*math.cos(angle1)*math.cos(angle2) \
            + z
        Y4 = (-lim3-float(self.probe_out.get()))*math.cos(angle1) \
            + (lim6+float(self.probe_out.get()))*math.sin(angle1) \
            + y
        X4 = (-lim1-float(self.probe_out.get()))*math.cos(angle2) \
            - (-lim3-float(self.probe_out.get()))*math.sin(angle1)*math.sin(angle2) \
            + (lim6+float(self.probe_out.get()))*math.cos(angle1)*math.sin(angle2) \
            + x

        # Create box parameters
        boxparam = u"""\n\t[SETTINGS.visiblebox]
\t# Coordinates of the vertices that define the visible 3D grid. Only four points are required to define the search space.
\n\t[SETTINGS.visiblebox.p1]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.visiblebox.p2]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.visiblebox.p3]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.visiblebox.p4]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.internalbox]
\t# Coordinates of the internal 3D grid. Used for calculations.
\n\t[SETTINGS.internalbox.p1]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.internalbox.p2]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.internalbox.p3]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
\n\t[SETTINGS.internalbox.p4]
\tx = {:.2f}
\ty = {:.2f}
\tz = {:.2f}
""".format(round(bX1, 2), round(bY1, 2), round(bZ1, 2),
           round(bX2, 2), round(bY2, 2), round(bZ2, 2),
           round(bX3, 2), round(bY3, 2), round(bZ3, 2),
           round(bX4, 2), round(bY4, 2), round(bZ4, 2),
           round(X1, 2), round(Y1, 2), round(Z1, 2),
           round(X2, 2), round(Y2, 2), round(Z2, 2),
           round(X3, 2), round(Y3, 2), round(Z3, 2),
           round(X4, 2), round(Y4, 2), round(Z4, 2))

        # Open parameters.toml
        f = open("parameters.toml", "w")
        # Write parameters in file
        f.write(param)
        f.write(boxparam)
        # Close parameters.toml
        f.close()

        return True

    def get_number_of_cavities(self):
        results = toml.load('{}/KV_Files/{}/{}.KVFinder.results.toml'.format(self.output.get(),
                                                                             self.base_name.getvalue(),
                                                                             self.base_name.getvalue()))
        return len(results['RESULTS']['VOLUME'].keys())

    def parKVFinder(self):
        """
        Run parKVFinder with parameters set in PyMOL parKVFinder Tools.
        :return: Execution success. True for success, False otherwise.
        :rtype: bool
        """
        # Running KVFinder
        sys.stdout.write(
            "\nRunning parKVFinder for: {}/KV_Files/{}.KVFinder.input.pdb\n".format(self.output.get(),
                                                                                  self.base_name.getvalue()))
        start_time = time.time()
        subprocess.call(self.execKVFinder.get().replace(' ', '\\ '),
                        stdout=subprocess.PIPE)
        number_of_cavities = self.get_number_of_cavities()
        sys.stdout.write("Number of cavities detected: {}!\n".format(number_of_cavities))
        elapsed_time = time.time() - start_time
        sys.stdout.write("done!\n")
        sys.stdout.write("Elapsed time: {:.2f} seconds\n\n".format(elapsed_time))

        # Copy parameters file in results folder and remove parameters file
        subprocess.call(['cp',
                         'parameters.toml',
                         '{}/KV_Files/parameters_{}.toml'.format(self.output.get(),
                                                                 self.base_name.getvalue())
                         ])
        subprocess.call(['rm',
                         'parameters.toml'])

        # Focus on Results Visualization
        self.notebook.selectpage('Results Visualization')
        self.notebook.tab('Results Visualization').focus_set()

        # Update Results File Entry
        self.results_file.configure(state='normal')
        self.results_file.delete(0, len(self.results_file.get()))
        self.results_file.insert(0,
                                 '{}/KV_Files/{}/{}.KVFinder.results.toml'.format(self.output.get(),
                                                                                  self.base_name.getvalue(),
                                                                                  self.base_name.getvalue()))
        self.results_file.configure(state='readonly')

        # Check if KVFinder find cavities
        if number_of_cavities > 0:
            # Save last objects used
            self.last_object.set('{}.KVFinder.output'.format(self.base_name.getvalue()))
            self.last_input.set('{}'.format(self.pdb_file.get(self.pdb_file.curselection())))
            if self.ligand_mode.get():
                self.last_ligand.set('{}'.format(self.ligand_file.get(self.ligand_file.curselection())))
            else:
                self.last_ligand.set('')

            # Load results
            self.load_results(0)

            return True
        elif number_of_cavities == 0:
            tkMessageBox.showwarning("Warning!", "No cavities found!", parent=self.parent)
            self.load_results(0)
            return True
        else:
            tkMessageBox.showerror("Error!", "Occurred an error in parKVFinder!", parent=self.parent)
            return False

    def show_grid(self):
        """
        Get minimum and maximum coordinates of the parKVFinder 3D-grid based on the set parameters.
        :return: Call draw_grid function with minimum and maximum coordinates or return Error.
        """
        global x, y, z

        if len(self.pdb_file.curselection()) > 0:
            # Get minimum and maximum dimensions of target PDB
            pdb = self.pdb_file.get(self.pdb_file.curselection())
            ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(pdb)

            # Get grid spacing method
            # Step size method (Grid spacing direct method)
            if self.grid_method.get():
                h = float(self.step_size.get())
            # Resolution method (Grid spacing indirect method)
            else:
                if self.resolution.getcurselection() == 'Low':
                    h = 0.6
                elif self.resolution.getcurselection() == 'Medium':
                    h = 0.5
                elif self.resolution.getcurselection() == 'High':
                    h = 0.25
            # Round dimensions
            if h == 0.0:
                tkMessageBox.showerror("Error", "Cannot draw grid! Step size set to 0.0!", parent=self.parent)
                return
            else:
                minX = round(minX - (minX % h), 1)
                minY = round(minY - (minY % h), 1)
                minZ = round(minZ - (minZ % h), 1)
                maxX = round(maxX - (maxX % h) + h, 1)
                maxY = round(maxY - (maxY % h) + h, 1)
                maxZ = round(maxZ - (maxZ % h) + h, 1)
            # Add probe_out
            probe_out = float(self.probe_out.get())
            probe_out = round(probe_out - round(probe_out, 4) % round(h, 4), 1)
            minX -= probe_out
            minY -= probe_out
            minZ -= probe_out
            maxX += probe_out
            maxY += probe_out
            maxZ += probe_out

            # Get center of each dimension (x, y, z)
            x = (minX + maxX)/2
            y = (minY + maxY)/2
            z = (minZ + maxZ)/2

            # Draw Grid
            self.draw_grid(minX, maxX, minY, maxY, minZ, maxZ)
        else:
            tkMessageBox.showerror("Error", "Load a PDB file!", parent=self.parent)
            return

    def draw_grid(self, minX, maxX, minY, maxY, minZ, maxZ):
        """
        Draw Grid in PyMOL.
        :param minX: minimum X coordinate.
        :param maxX: maximum X coordinate.
        :param minY: minimum Y coordinate.
        :param maxY: maximum Y coordinate.
        :param minZ: minimum Z coordinate.
        :param maxZ: maximum Z coordinate.
        :return: grid object in PyMOL.
        """
        angle1 = 0
        angle2 = 0
        lim1 = x - minX
        lim2 = maxX - x
        lim3 = y - minY
        lim4 = maxY - y
        lim5 = z - minZ
        lim6 = maxZ - z

        # Get positions of box vertices
        z1 = lim1*math.sin(float(angle2)) \
            + lim3*math.sin(float(angle1))*math.cos(float(angle2)) \
            - lim5*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y1 = -lim3*math.cos(float(angle1)) \
            + (-lim5)*math.sin(float(angle1)) \
            + y
        x1 = -lim1*math.cos(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z2 = (-lim2)*math.sin(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.cos(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y2 = (-lim3)*math.cos(float(angle1)) \
            + (-lim5)*math.sin(float(angle1)) \
            + y
        x2 = lim2*math.cos(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z3 = -(-lim1)*math.sin(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.cos(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y3 = lim4*math.cos(float(angle1)) \
            + (-lim5)*math.sin(float(angle1)) \
            + y
        x3 = (-lim1)*math.cos(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.sin(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z4 = -(-lim1)*math.sin(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.cos(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y4 = (-lim3)*math.cos(float(angle1)) \
            + lim6*math.sin(float(angle1)) \
            + y
        x4 = (-lim1)*math.cos(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z5 = (-lim2)*math.sin(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.cos(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y5 = lim4*math.cos(float(angle1)) \
            + (-lim5)*math.sin(float(angle1)) \
            + y
        x5 = lim2*math.cos(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.sin(float(angle2)) \
            + (-lim5)*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z6 = (-lim2)*math.sin(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.cos(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y6 = (-lim3)*math.cos(float(angle1)) \
            + lim6*math.sin(float(angle1)) \
            + y
        x6 = lim2*math.cos(float(angle2)) \
            - (-lim3)*math.sin(float(angle1))*math.sin(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z7 = -(-lim1)*math.sin(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.cos(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y7 = lim4*math.cos(float(angle1)) \
            + lim6*math.sin(float(angle1)) \
            + y
        x7 = (-lim1)*math.cos(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.sin(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x
        z8 = (-lim2)*math.sin(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.cos(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.cos(float(angle2)) \
            + z
        y8 = lim4*math.cos(float(angle1)) \
            + lim6*math.sin(float(angle1)) \
            + y
        x8 = lim2*math.cos(float(angle2)) \
            - lim4*math.sin(float(angle1))*math.sin(float(angle2)) \
            + lim6*math.cos(float(angle1))*math.sin(float(angle2)) \
            + x

        # Create box object
        if "grid" in cmd.get_names("objects"):
            cmd.delete("grid")

        # Create vertices
        cmd.pseudoatom("grid", name="v2", pos=[x2, y2, z2], color="white")
        cmd.pseudoatom("grid", name="v3", pos=[x3, y3, z3], color="white")
        cmd.pseudoatom("grid", name="v4", pos=[x4, y4, z4], color="white")
        cmd.pseudoatom("grid", name="v5", pos=[x5, y5, z5], color="white")
        cmd.pseudoatom("grid", name="v6", pos=[x6, y6, z6], color="white")
        cmd.pseudoatom("grid", name="v7", pos=[x7, y7, z7], color="white")
        cmd.pseudoatom("grid", name="v8", pos=[x8, y8, z8], color="white")

        # Connect vertices
        cmd.select("vertices", "(name v3,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v5,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v2,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v6)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v4,v7)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v3,v5)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v6,v8)")
        cmd.bond("vertices", "vertices")
        cmd.select("vertices", "(name v7,v8)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("grid", name="v1x", pos=[x1, y1, z1], color='white')
        cmd.pseudoatom("grid", name="v2x", pos=[x2, y2, z2], color='white')
        cmd.select("vertices", "(name v1x,v2x)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("grid", name="v1y", pos=[x1, y1, z1], color='white')
        cmd.pseudoatom("grid", name="v3y", pos=[x3, y3, z3], color='white')
        cmd.select("vertices", "(name v1y,v3y)")
        cmd.bond("vertices", "vertices")
        cmd.pseudoatom("grid", name="v4z", pos=[x4, y4, z4], color='white')
        cmd.pseudoatom("grid", name="v1z", pos=[x1, y1, z1], color='white')
        cmd.select("vertices", "(name v1z,v4z)")
        cmd.bond("vertices", "vertices")
        cmd.delete("vertices")

    def restore_values(self):
        """
        Restore all parameters and results to their default values.
        """
        clean_results = tkMessageBox.askyesno("Restore Values", "Also restore Results Visualization tab?")
        if clean_results:
            ## Remove cavities, residues and last_object objects
            cmd.delete("cavities")
            cmd.delete("residues")
            cmd.delete(self.last_object.get())
            self.last_object.set("")
            if self.loaded_results.get():
                cmd.delete(self.last_input.get())
                cmd.delete(self.last_ligand.get())
                self.last_input.set("")
                self.last_ligand.set("")
            cmd.frame(1)
            self.clean_results()

            # Remove Results File Entry
            self.results_file.configure(state='normal')
            self.results_file.delete(0, len(self.results_file.get()))
            self.results_file.configure(state='readonly')

        #### Restore parameters
        cmd.delete("grid")
        self.surface.setvalue(self.defaults['surface_representation'])
        self.cavity_representation.setvalue(self.defaults['cavity_representation'])
        self.base_name.delete(0, len(self.base_name.get()))
        self.base_name.insert(0, self.defaults['base_name'])
        self.probe_in.setentry(self.defaults['probe_in'])
        self.probe_out.setentry(self.defaults['probe_out'])
        self.volume_cutoff.setentry(self.defaults['volume_cutoff'])
        self.removal_distance.setentry(self.defaults['removal_distance'])
        self.output.configure(state='normal')
        self.output.delete(0, len(self.output.get()))
        self.output.insert(0, self.wd.get())
        self.output.configure(state='readonly')
        self.grid_method.set(self.defaults['grid_method'])
        self.step_size.setentry(self.defaults['step_size'])
        self.resolution.setvalue(self.defaults['resolution'])
        self.search_procedure.setvalue(self.defaults['search_procedure'])
        self.update_box_adjustment_frame(self.defaults['search_procedure'])
        self.ligand_mode.set(self.defaults['ligand_mode'])
        self.ligand_cutoff.setentry(self.defaults['ligand_cutoff'])
        self.check_ligand_mode()
        cmd.frame(1)

    def quit(self):
        """
        Destroy PyMOL parKVFinder Tools Pmw.Dialog.
        """
        self.dialog.destroy()

    #### (6.7) Function to handle execution process in the Pmw.Dialog ##################################################
    def execute(self, result, refocus=True):
        """
        Handle execution process of PyMOL parKVFinder Tools Pmw.Dialog.
        :param result: Button name.
        :param refocus: Focus Pmw.Dialog.
        :return: Output of each button called.
        """
        if result == 'Run parKVFinder':
            # Set working directory for path defined in Output Directory Entry
            os.chdir(self.output.get())
            # Save parameters file
            good = self.create_parKVFinder_input()
            if not good:
                tkMessageBox.showerror("Error",
                                       "Cannot generate parameters file. Verify chosen parameters!",
                                       parent=self.parent)
                return False
            # Run parKVFinder using created parameters file
            self.parKVFinder()
        elif result == 'Show Grid':
            self.show_grid()
        elif result == 'Restore Default Values':
            self.restore_values()
        elif result == 'Save Parameters':
            # Set working directory for path defined in Output Directory Entry
            os.chdir(self.output.get())
            # Save parameters file
            good = self.create_parKVFinder_input()
            if not good:
                tkMessageBox.showerror("Error",
                                       "Cannot generate parameters file. Verify chosen parameters!",
                                       parent=self.parent)
                return False
        else:
            self.quit()

    #### (6.8) Show Dialog App Function ################################################################################
    def show_app_modal(self):
        """
        Display PyMOL parKVFinder Tools Pmw.Dialog.
        :return: Display continuosly Pmw.Dialog of PyMOL parKVFinder Tools.
        """
        self.dialog.show()
