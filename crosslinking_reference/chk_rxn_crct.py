"""
This module houses the CheckReactionCrosslink class which performs the correcting routine as described in Freger_2014_JMemSci.

Imports the following polymerizeit modules:
    - polymerizeit.reaction.custom.correcting_routine.reaction_crietria_functions
    - polymerizeit.reaction.reaction_criteria.BondPiercing
    - polymerizeit.reaction.core.merge_gro_tables
    - polymerizeit.reaction.additional_features_to_core.insrt_mol
    - polymerizeit.reaction.custom.correcting_routine.calc_dist

:Author(s): Siva Dasetty, Salman Bin Kashif


:Version History:
    **2021-01-27**: Created the module
    **2022-12-XX**: Repurposed the module for the new architecture of polymerizeit
    **2023-12-XX**: Updated to (another) new architecture of polymerizeit core
"""

# standard python modules
import sys
import sqlite3
import os
import re
import math
import importlib

import polymerizeit as pi
from polymerizeit.reaction.custom.correcting_routine import reaction_crietria_functions as chk
from polymerizeit.reaction.reaction_criteria.BondPiercing import BondPiercing

class CheckReactionCrosslink(pi.ProtocolBase):
    """
    
    Implements the correcting routine as described in Freger_2014_JMemSci. May be repurposed later as standard crosskinking routine

    :param pi.ProtocolBase: Inherits from the ProtocolBase class
    :type pi.ProtocolBase: class
        
    """
    def __init__(self,**kwargs):
        """Standard module to check crosslinking reactions
        
        This is a standard module to check crosslinking reactions.
        
        """
        super().__init__(**kwargs)
        self.__dict__.update(kwargs)
        self.mon_A_name_r=[]
        self.mon_B_name_r=[]
        self.clust_A_r=[]
        self.clust_B_r=[]
        self.tri_A_r=[]
        self.tri_B_r=[]
        self.blocked_clust_pairs={}
        self.rx_corct=getattr(importlib.import_module("polymerizeit.reaction.custom.correcting_routine.rxn_crct"),"ReactionCrosslink")(preprocessed_inputs=self.preprocessed_inputs)
        self.chk=importlib.import_module("polymerizeit.reaction.custom.correcting_routine.calc_dist")
        self.cg=getattr(importlib.import_module("polymerizeit.reaction.core.merge_gro_tables"),"UpdateStructure")(preprocessed_inputs=self.preprocessed_inputs)
        self.im=getattr(importlib.import_module("polymerizeit.reaction.additional_features_to_core.insrt_mol"),"InsertMolecule")(preprocessed_inputs=self.preprocessed_inputs)
        self.grofile=self.ifgro
    
    def calculate_tau(self): 
        """
        Measures the deviation of the system from the target atomic composition
        
        Parameters
        ----------
        None
        
        Returns
        -------
        avg_tau: float
            Value of tau for the input system
        """
        
        #print("Module chk_rxn_crct.py. Expected wet-lab composition:", self.exp_vals)

        # calculate tau for each cluster in the system.
        sqldb = self.main_database
        conn = sqlite3.connect(sqldb)
        c = conn.cursor()
        
        # master table that contains information of all chains and monomers 
        # fetch table name
        tabmol = c.execute("""SELECT name FROM sqlite_master WHERE type = ? AND name like ?;""", ('table', '%system%mol%', )).fetchone()[0]

        # master table that contains information about the number of carbon, oxygen, nitrogen, chlorine, and hydrogen atoms in
        # molecule in the system table 
        # fetch table name
        tabatm = c.execute("""SELECT name FROM sqlite_master WHERE type = ? AND name like ?;""", ('table', '%system%atm%', )).fetchone()[0]

        # here clusters indicate both monomers and polymer chains; list contains names of clusters in the system 
        # monomers are excluded while reading the list (see below)
        clusters = c.execute('SELECT molname FROM {tab} WHERE nmol != ?'.format(tab=tabmol), (0,)).fetchall()

        #if self.inputs_file.split(".")[-1]=="txt":
        self.exp_vals=self.exp_vals.split(',')
        
        Ac_exp = float(self.exp_vals[0])
        Ao_exp = float(self.exp_vals[1])
        An_exp = float(self.exp_vals[2])
        Acl_exp = float(self.exp_vals[3])
        Aon_exp = float(self.exp_vals[4])

        #print("Module chk_rxn_crct. Function calculate_tau. Varaiable Ac_exp", Ac_exp)
    
        # read tabatm, fetch num_c, num_o, num_n, num_h, and num_cl
        tau = 0
        count_clusters = 0
        #print(clusters)
        for cluster in clusters:
            #print("Module check_rxn_crct.py. Checking composition for:", cluster[0])
            # exlude monomers;
            if cluster[0] == self.mon_A_name or cluster[0] == self.mon_B_name:
                continue
            # get num_atm values from tabatm, which keeps track of number of each atom in a cluster/monomer
            num_c = int(c.execute('SELECT num_c FROM {tab} WHERE molname == ?'.format(tab=tabatm), (cluster[0],)).fetchone()[0])
            num_o = int(c.execute('SELECT num_o FROM {tab} WHERE molname == ?'.format(tab=tabatm), (cluster[0],)).fetchone()[0])
            num_n = int(c.execute('SELECT num_n FROM {tab} WHERE molname == ?'.format(tab=tabatm), (cluster[0],)).fetchone()[0])
            num_cl = int(c.execute('SELECT num_cl FROM {tab} WHERE molname == ?'.format(tab=tabatm), (cluster[0],)).fetchone()[0])
            count_clusters += 1
            # total number
            tot_atms = num_c + num_o + num_n + num_cl

            # calculate Ac, Ao, An, Acl, Aon for each chain in the system
            Ac = float(num_c/tot_atms)
            Ao = float((num_o + num_cl)/tot_atms)
            An = float(num_n/tot_atms)
            Acl = float(num_cl/tot_atms)
            Aon = float(Ao/An)

            # calculate sum of abs of tau's of all chains in the system;
            tau += float(format(abs(5.000 - Ac/Ac_exp + Ao/Ao_exp +  An/An_exp +  Acl/Acl_exp +  Aon/Aon_exp), '8.3f'))

        #print("Module chk_rxn_crct.py. Value of tau is:", tau)
        # average of tau
        avg_tau = tau/count_clusters

        return avg_tau


    # check variation of tau with cycle every 5 cycles
    # if variation is monotonic, return 1 --> enable correcting routine
    def check_tau(self,
                  tau_file):
        """
        Check the variation of tau in each cycle every 5 cycles

        Parameters
        ----------
        tau_file: String
            File recording the tau value from previous cycles
        
        Returns
        -------
        corct_routine: int
            Binary output. 1 if correcting is to be enabled, otherwise 0.

        """
        # if everything went right; there should be five lines in tau.txt
        try:
            f = open(tau_file, 'r')
        except IOError as e:
            print("Unable to open " + tau_file + ". Please check the file.")

        # check whether the value in current line is higher than previous line
        tau_val_prev = -1
        for line in f:
            tau_val_after = float(line) - tau_val_prev

            if tau_val_after > 0:
                tau_val_prev =  float(line)
            else:
                break


        corct_routine = 0

        # if the value is still positive; it means monotonic increase
        if tau_val_after > 0:
            corct_routine = 1

        f.close()

        return corct_routine
    def check_tau_dict(self,
                       tau_dict,
                       iteration):
        """
        Check the variation of tau in each cycle every 5 cycles
        
        Parameters
        ----------
        tau_dict: dict
            Dictionary containing the tau values from previous cycles
        iteration: int
            Current iteration number
        Returns
        -------
        corct_routine: int
            Binary output. 1 if correcting is to be enabled, otherwise 0.
        """
        #Obatain the tau values for the previous five iterations
        tau_list=[]
        for i in range(iteration-4,iteration+1):
            tau_list.append(tau_dict[i])
        #Check if the tau values are monotonously increasing
        if all(round(tau_list[i],5) < round(tau_list[i+1],5) for i in range(len(tau_list)-1)):
            chk_tau=1
            
        else:
            chk_tau=0
            
        return chk_tau
                 


    # Correcting routine algorithm: Control COCL concentration in the system by cross-linking.
    #     Find COCL sites in the polymer chains that are within 2.5 Angstroms, add MPD monomer between the sites
    #     to form two bonds simultaneously with the two COCL sites. If there is already a MPD monomer, perform reaction
    #     with that monomer. Here, if the two sites are part of the same chain, then the two sites should be separated 
    #     by at least n (default = 10) monomers. Add the MPD monomer after finding the target orientation. 
    #     In case, there are no COCL sites within 2.5 Angstroms, find sites within 5 Angstroms. Here again, if the sites are
    #     part of the same chain, they should be separated by at least 10 monomers. Then, add a MPD-TMC-MPD trimer (if there is 
    #     no trimer in that location) between the sites. Create bonds between the COCL sites and MPD sites of the trimer. 
    #     All cutoff values are user defined.

