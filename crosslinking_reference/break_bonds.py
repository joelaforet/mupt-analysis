"""Generate the sqlite tables from top and gro files

This module houses the GenMSTabs class. This class deals with creating the
sqlite databases from the input files, which are needed later by ReactionStandard 
class and some custom classes to implement the topology change after the reaction.

Authors
-------
Salman Bin Kashif, Siva Dasetty, Ryan DeFever, Sapna Sarupria

Example
-------
Not needed


Versions
--------
Created: 01 Oct 2017
Edited: 22 Jan 2018; Updated comments
Edited: 15 Jun 2018; Updated gen_dih to include improper terms as well. Check for (#update_20180715)
Edited: 17 March 2020; Modified the filenames to work with generalized code
Edited: May 1, 2021; Modified to work with multiple repeat units
Edited: May 2022: Modified to work with new (automated) pre-processing scheme
"""



# standard python modules
import sys
import re
import sqlite3
import math
import os
import inspect
from sqlite3 import Error


class GenSystemTables(object):
    """Creates the sqlite tables from the top and gro files at each iteration

    A polymerization procotocol typically involves multiple iterations of finding the
    new bonds, updating the topology as per the new bonds, and then equilibrating the system.
    This class is called at every new iteration to read the topology information from last iteration
    after the equilibraion steps. 
    """
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)

    def from_gromacs(self,function="main"):
        """
        from_gromacs _summary_

        _extended_summary_

        Parameters
        ----------
        function : str, optional
            _description_, by default "main"

        Returns
        -------
        _type_
            _description_
        """        
        print("Reading the gro and top files into database ...")
        print(self.ifgro,self.iftop,self.main_database)
        conn = sqlite3.connect(self.main_database)
        c=conn.cursor()
        x = c.execute("""SELECT name FROM sqlite_master WHERE type = ?;""", ('table', )).fetchall() 
    
        table_names = map(lambda tup: tup[0], x)
        
        for name in table_names:
            c.execute("DROP TABLE IF EXISTS %s" % name)
        conn.commit()
        conn.close()
        #self.create_gro_table_from_gromacs(grofile=self.ifgro)
        #self.create_top_table_from_gromacs(topfile=self.iftop,function=function)
        #Delete database if it exists
        os.remove(self.main_database)
        self.create_itptab_from_gromacs(self.itpfile)
        print("Database created...")
        
        base = os.path.basename(self.itpfile)
        keyfname = os.path.splitext(base)[0]
        keyfname=keyfname.split('_')[0]
        
        keyfname=keyfname.split('-')[0]
        target_dc=self.target_dc
        
        self.break_bonds(keyfname,target_dc)
        self.write_itpfile(self.itpfile,self.main_database,keyfname,target_dc)
        return None



    # create tables for atm, bnd, ang, and dih using .itp file. 
    def create_itptab_from_gromacs(self,itpfile):
        """_summary_

        Parameters
        ----------
        itpfile : _type_
            _description_
        """        
    
        # Open .itp file
        try:
            f = open(itpfile, 'r')
        except IOError as e:
            print("Unable to open " + itpfile + ". Please check the file.")

        base = os.path.basename(itpfile)
        keyfname = os.path.splitext(base)[0]
        keyfname=keyfname.split('_')[0]
        
        keyfname=keyfname.split('-')[0]

        #print("Module gen_molstabs.py. Function create_itptab_from_gromacs. Variable keyfname.", keyfname)

    	# connect to sqlite database and assign table names
        sqldb = self.main_database
        conn = sqlite3.connect(sqldb)
        c = conn.cursor()
        
        #clear the tables if they exist
        c.execute('drop table if exists {tab}'.format(tab=keyfname+"_pi_master_atm"))
        c.execute('drop table if exists {tab}'.format(tab=keyfname+"_pi_master_bnd"))
        c.execute('drop table if exists {tab}'.format(tab=keyfname+"_pi_master_ang"))
        c.execute('drop table if exists {tab}'.format(tab=keyfname+"_pi_master_dih"))

        tabatm = keyfname + '_pi_master_atm' 
        tabbnd = keyfname + '_pi_master_bnd' 
        tabang = keyfname + '_pi_master_ang' 
        tabdih = keyfname + '_pi_master_dih' 

        # create tables
        c.execute('CREATE TABLE IF NOT EXISTS {tab} (ai text, name text, resid int, resnm text, atname text, cgnr int, charge double, mass real, uniq_atname text, mnm_type text)'.format(tab=tabatm))
        c.execute('CREATE TABLE IF NOT EXISTS {tab} (ai text, aj text, funct int, c0 real, bondc float)'.format(tab=tabbnd))
        c.execute('CREATE TABLE IF NOT EXISTS {tab} (ai text, aj text, ak text, funct int, c0 real, c1 real)'.format(tab=tabang))
        c.execute('CREATE TABLE IF NOT EXISTS {tab} (ai text, aj text, ak text, al text, funct int, c0 real, c1 real, mult int)'.format(tab=tabdih)) 


        # modify here if pairs section is also read
        # p*_pi.itp files may not have pairs section
  

        psec = "[ pairs ]"
        pflag = 3
        asec = "[ angles ]"
        aflag = 4


        # read into tables; pairs section ignored here -- use gen-pairs in .top file 
        flag = 0
        for line in f:
            # comments in .itp file start with ";" and ignore empty lines
            if not line.startswith(";") and line.strip():

                if line.strip() == "[ atoms ]":
                    flag = 1
                    continue

                elif flag == 1 and not line.strip() == "[ bonds ]":
                    linedata = line.strip().split()[:13]
                    atmvalues = (linedata[0], linedata[1], linedata[2], linedata[3], linedata[4], linedata[5], linedata[6], linedata[7], linedata[11], linedata[12])
                    c.execute('INSERT INTO {tab} VALUES(?,?,?,?,?,?,?,?,?,?)'.format(tab=tabatm), atmvalues)

                elif line.strip() == "[ bonds ]":
                    flag = 2
                    continue

                elif flag == 2 and not line.strip() == psec:
                    linedata = line.strip().split()
                    if len(linedata) == 5:
                        bndvalues = (linedata[0], linedata[1], linedata[2], linedata[3], linedata[4])
                    elif len(linedata) == 3:
                        bndvalues = (linedata[0], linedata[1], linedata[2], "", "")
                    else:
                        print("Bond list appears to be incomplete. Please check " + itpfile + ". \n")
                        exit()
                    c.execute('INSERT INTO {tab} VALUES(?,?,?,?,?)'.format(tab=tabbnd), bndvalues)

                elif line.strip() == psec:
                    flag = pflag
                    continue

                elif line.strip() == asec:
                    flag = aflag
                    continue

                elif flag == 4 and not line.strip() == "[ dihedrals ]":
                    linedata = line.strip().split()
                    if len(linedata) == 6:
                        angvalues = (linedata[0], linedata[1], linedata[2], linedata[3], linedata[4], linedata[5])
                    elif len(linedata) == 4:
                        angvalues = (linedata[0], linedata[1], linedata[2], linedata[3], "", "")
                    else:
                        print("Angle list appears to be incomplete. Please check " + itpfile + ". \n")
                        exit()
                    c.execute('INSERT INTO {tab} VALUES(?,?,?,?,?,?)'.format(tab=tabang), angvalues)

                elif line.strip() == "[ dihedrals ]":
                    flag  = 5
                    continue

                elif flag == 5 and not line.strip() == "[ system ]":
                    linedata = line.strip().split()
                    if len(linedata) == 8:
                        dihvalues = (linedata[0], linedata[1], linedata[2], linedata[3], linedata[4], linedata[5], linedata[6], linedata[7])
                    elif len(linedata) == 5:
                        dihvalues = (linedata[0], linedata[1], linedata[2], linedata[3], linedata[4], "", "", "")
                    else:
                        print("Dihedral list appears to be incomplete. Please check " + itpfile + ". \n")
                        exit()
                    c.execute('INSERT INTO {tab} VALUES(?,?,?,?,?,?,?,?)'.format(tab=tabdih), dihvalues)

                elif line.strip() == "[ system ]":
                    flag = 6
                    continue

        conn.commit()
        conn.close()
        f.close()
        
    def break_bonds(self,keyfname,target_dc=0.5):
        
        print("Breaking bonds to achieve target degree of crosslinking...")
        #Read in the database
        conn = sqlite3.connect(self.main_database)
        c = conn.cursor()
        
        # tables
        tabatm = keyfname + '_pi_master_atm'
        tabbnd = keyfname + '_pi_master_bnd'
        tabang = keyfname + '_pi_master_ang'
        tabdih = keyfname + '_pi_master_dih'
        
        import pandas as pd
        #Use sqlite qrey to  read in the atom indices of NB and CB
        # ai_n_type=pd.read_sql_query(f"select ai from {tabatm} where atname like 'NB%'",con=conn)
        # ai_c_type=pd.read_sql_query(f"select ai from {tabatm} where atname like 'CB%'",con=conn)
        
        ai_n_type=pd.read_sql_query(f"select * from {tabatm} where atname like 'NB%' and resnm like 'MPD2'",con=conn)
        ai_c_type=pd.read_sql_query(f"select ai from {tabatm} where atname like 'CB%' and (resnm like 'TMC1' or resnm like 'TMC2' or resnm like 'TMC3')",con=conn)
        
        ai_n_type_all=pd.read_sql_query(f"select ai from {tabatm} where atname like 'NB%'",con=conn)
        
        
        ai_c_type_list=ai_c_type['ai'].tolist()
        
        # Apply the filter with groupby and apply
        filtered_n_type = ai_n_type.groupby(['resid', 'resnm'], group_keys=False).apply(self.random_filter).reset_index(drop=True)

        
        
        breakpoint()
        
        
        #assert that both reacted atoms are same in number
        #assert len(ai_n_type)==len(ai_c_type), "Number of atoms in the two types are not same"
        
        #Calculate total number of N and C atoms
        total_n_atoms=pd.read_sql_query(f"select count(ai) from {tabatm} where atname like 'N%'",con=conn).values[0][0]
        total_c_atoms=pd.read_sql_query(f"select count(ai) from {tabatm} where atname like 'C%' and not atname like 'CA%'",con=conn).values[0][0]
        
        #Calculate the degree of crosslinking -- n/(n+min((n or c)))
        dc=ai_n_type_all.shape[0]/(min(total_n_atoms,total_c_atoms))
        
        breakpoint()
        n_bonds_to_break=math.ceil((dc-target_dc)*min(total_n_atoms,total_c_atoms))
        
        #If number of bonds to break is less than filterd_n_type, then sample from filtered_n_type
        if n_bonds_to_break<filtered_n_type.shape[0]:
            import random
            # Sample one-time unique indices for NB
            sampled_n_indices = random.sample(filtered_n_type['ai'].tolist(), n_bonds_to_break)  # Adjust the sample size as needed
        else:
            breakpoint()
            #first tbe the possible bonds to break, thr pick the dfference from MPI
            sampled_n_indices=filtered_n_type['ai'].tolist()
            n_bonds_to_break=n_bonds_to_break-len(sampled_n_indices)
            ai_mpd1_n_type=pd.read_sql_query(f"select * from {tabatm} where atname like 'NB%' and resnm like 'MPD1'",con=conn)
            # Sample one-time unique indices for NB
            sampled_n_indices = random.sample(ai_mpd1_n_type['ai'].tolist(), n_bonds_to_break)  # Adjust the sample size as needed
            breakpoint()
        # Find bonds to break
        query_bonds = f"""
        SELECT ai, aj FROM {tabbnd}
        WHERE ai IN ({",".join(map(str, sampled_n_indices))}) 
        AND aj IN ({",".join(map(str, ai_c_type_list))})
        OR (ai IN ({",".join(map(str, ai_c_type_list))}) AND aj IN ({",".join(map(str, sampled_n_indices))}))
        """
        bond_pairs = pd.read_sql_query(query_bonds, con=conn)

        deleted_bonds = []
        deleted_angles = []
        deleted_dihedrals = []
        
        
        if bond_pairs.empty:
            print("No bond pairs found to break.")
        else:
            for _, row in bond_pairs.iterrows():
                print(f"Processing bond pair: {row['ai']} - {row['aj']}")
                ai, aj = row['ai'], row['aj']
                deleted_bonds.append((ai, aj))

                # Query and delete angles
                query_angles = f"SELECT ai, aj, ak FROM {tabang} WHERE (ai = {ai} AND aj = {aj}) OR (ai = {ai} AND ak = {aj}) OR (aj = {ai} AND ak = {aj}) OR (ai = {aj} AND aj = {ai}) OR (ai = {aj} AND ak = {ai}) OR (aj = {aj} AND ak = {ai})"
                angles_to_delete = pd.read_sql_query(query_angles, con=conn)
                for _, angle_row in angles_to_delete.iterrows():
                    deleted_angles.append((angle_row['ai'], angle_row['aj'], angle_row['ak']))
                    c.execute(f"DELETE FROM {tabang} WHERE ai = ? AND aj = ? AND ak = ?", (angle_row['ai'], angle_row['aj'], angle_row['ak']))

                # Query and delete dihedrals
                query_dihedrals = f"SELECT ai, aj, ak, al FROM {tabdih} WHERE (ai = {ai} AND aj = {aj}) OR (ai = {ai} AND ak = {aj}) OR (ai = {ai} AND al = {aj}) OR (aj = {ai} AND ak = {aj}) OR (aj = {ai} AND al = {aj}) OR (ak = {ai} AND al = {aj}) OR (ai = {aj} AND aj = {ai}) OR (ai = {aj} AND ak = {ai}) OR (ai = {aj} AND al = {ai}) OR (aj = {aj} AND ak = {ai}) OR (aj = {aj} AND al = {ai}) OR (ak = {aj} AND al = {ai})"
                dihedrals_to_delete = pd.read_sql_query(query_dihedrals, con=conn)
                for _, dihedral_row in dihedrals_to_delete.iterrows():
                    deleted_dihedrals.append((dihedral_row['ai'], dihedral_row['aj'], dihedral_row['ak'], dihedral_row['al']))
                    c.execute(f"DELETE FROM {tabdih} WHERE ai = ? AND aj = ? AND ak = ? AND al = ?", (dihedral_row['ai'], dihedral_row['aj'], dihedral_row['ak'], dihedral_row['al']))

                # Delete the bond last to maintain integrity during angle and dihedral deletion
                c.execute(f"DELETE FROM {tabbnd} WHERE (ai = ? AND aj = ?) OR (ai = ? AND aj = ?)", (ai, aj, aj, ai))

            conn.commit()
            print(f"Bonds deleted: {deleted_bonds}")
            print(f"Angles deleted: {deleted_angles}")
            print(f"Dihedrals deleted: {deleted_dihedrals}")
            
            #Save to txt file
            with open(f"{keyfname}_deleted_bonds_{target_dc}.txt", "w") as f:
                for bond in deleted_bonds:
                    f.write(f"{bond}\n")
            with open(f"{keyfname}_deleted_angles.txt_{target_dc}", "w") as f:
                for angle in deleted_angles:
                    f.write(f"{angle}\n")
            with open(f"{keyfname}_deleted_dihedrals_{target_dc}.txt", "w") as f:
                for dihedral in deleted_dihedrals:
                    f.write(f"{dihedral}\n")

        conn.close()
        
        
    # Custom function to filter records randomly based on resnm
    def random_filter(self,group):
        if group['resnm'].iloc[0] in ['MPD2', 'TMC2']:
            # Randomly return one record for MPD2 or TMC2
            return group.sample(n=1)
        elif group['resnm'].iloc[0] == 'TMC3':
            # Randomly return up to two records for TMC3, checking if the group size allows for it
            return group.sample(n=min(2, len(group)))
        else:
            # Default behavior for other resnm values, if any
            return group
        
        
    def write_itpfile(self,itpfile,database,molecule_name,target_dc):
        """
        write_itpfile _summary_

        Parameters
        ----------
        itpfile_path : _type_
            _description_
        """        

        filepath = os.path.dirname(itpfile)
        filename = os.path.join(filepath, f"{molecule_name}_dc_{target_dc}_.itp")
        
        print("Writing update itp file to ",filename)

            # Continue with the rest of your code...
        f = open(filename, 'w')

        sqldb = database
        conn = sqlite3.connect(sqldb)
        c = conn.cursor()
        
        
        import getpass
        import socket
        from datetime import datetime

        username=getpass.getuser()
        hostname=socket.gethostname()
        system_info=os.uname()
        where=os.getcwd()
        f.write("; Generated by break_bonds.py\n")
        f.write("; User name:" + (str) (username) + "\n")
        f.write("; Host:"+(str) (hostname) +"\n")
        f.write("; System info:"+(str) (system_info) + "\n")
        f.write("; Current directory: "+ (str) (where)+"\n")
        f.write("; Created: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n\n")


        f.write("[ moleculetype ]" + "\n")
        f.write("; Name \t \t \t nrexcl" + "\n")
        f.write("Polyamide" + "\t \t \t" + str(3) + "\n \n")
    

        # print net_charge; warn user if system is not neutral
        tabatm = molecule_name.lower() + "_pi_master_atm"
        net_charge = c.execute('SELECT charge FROM {tab}'.format(tab = tabatm)).fetchall()

        # tot_charge
        tot_charge = 0
        for row in net_charge:
            tot_charge += float('{:.15f}'.format(row[0]))

        charge = 0.0 
        for sec in ["atm", "bnd", "ang", "dih"]:
           tabnm = molecule_name.lower() + "_pi_master_" + sec

           if sec == "atm":
               f.write("[ atoms ]" + "\n")
               f.write(";   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB    massB  unique_atoms mnm_type \n")
               f.write("; residue    1 DIM rtp DIM q  " + str(format(tot_charge, '8.6f')) + "\n")
               vals = c.execute('SELECT * from {tab}'.format(tab = tabnm)).fetchall()

           elif sec == "bnd":
               f.write("[ bonds ]" + "\n")
               f.write(";    ai     aj funct         c0         c1         c2         c3\n")
               vals = c.execute('SELECT * from {tab} ORDER BY funct'.format(tab = tabnm)).fetchall()

           elif sec == "prs":
               f.write("[ pairs ]" + "\n \n")
        #       vals=c.execute('SELECT * from {tab} ORDER by funct'.format(tab=tabnm)).fetchall()

           elif sec == "ang":
               f.write("[ angles ]" + "\n")
               f.write(";    ai     aj     ak funct         c0         c1         c2         c3\n")
               vals = c.execute('SELECT * from {tab} ORDER BY funct'.format(tab = tabnm)).fetchall()

           elif sec == "dih":
               f.write("[ dihedrals ]" + "\n")
               f.write(";    ai     aj     ak     al funct         c0         c1         c2         c3         c4         c5\n")
               vals = c.execute('SELECT * from {tab} ORDER by funct'.format(tab = tabnm)).fetchall()

           for row in  vals:
               if sec == "atm":
                   charge += float('{:.15f}'.format(row[6]))
                   f.write("\t" + row[0] + "\t" +  row[1] + "\t")
                   f.write(str(row[2]) + "\t" + row[3] + "\t") 
                   f.write(row[4] + "\t" + str(row[5]) + "\t")
                   f.write(str(format(row[6], '9.15f')) + "\t" + str(row[7]) + "\t")
                   #f.write("; "+ "qtot" + "\t" + str(format(charge, '9.6f')) + "\t" + format(row[11], '10s') + "\t" + "\t" + format(row[12], '5s') + "\t")
               elif sec == "prs":
                   continue
               else:
                   #> -- right alignment; 8 number of characters
                   # values are read from ffbonded.itp file
                   f.write(('{!s:>8}   '*len(row)).format(*row))

               f.write("\n")
           f.write("\n")

        f.write("\n \n \n")
        f.close()


    
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="Reads the gro and top files into database")
    parser.add_argument('--ifgro', type=str, help='Input gro file')
    parser.add_argument('--iftop', type=str, help='Input top file')
    parser.add_argument('--itpfile', type=str, help='Membrane itp file')
    parser.add_argument('--main_database', type=str, help='Main database file')
    parser.add_argument('--target_dc', type=float, help='Target degree of crosslinking', default=0.7)
    
    args = parser.parse_args()
    
    read_obj = GenSystemTables(ifgro=args.ifgro, 
                               iftop=args.iftop,
                               itpfile=args.itpfile, 
                               main_database=args.main_database,
                               target_dc=args.target_dc)
    read_obj.from_gromacs()
