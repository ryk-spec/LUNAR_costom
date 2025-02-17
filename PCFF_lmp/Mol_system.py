from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
import os,subprocess

packmol_path = "/home/ryo/projects/packmol/packmol"
class Molecule:
    def __init__(self, mol_obj : Chem.Mol, name : str):
        self.name = name
        self.mol_obj = mol_obj
        self.gram = Descriptors.MolWt(mol_obj)/6.022e23

    def load_pdb(self):

        num_conformers = 1
        conformer_ids = AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=num_conformers)
        # PDBブロックとして各コンフォーマーを保存
        pdb_block = Chem.MolToPDBBlock(self.mol_obj, confId=0)
        pdb_block = pdb_block.replace("BR","Br")
        pdb_block = pdb_block.replace("SI","Si")
        pdb_name = f"{self.name}.pdb"
        with open(pdb_name, "w") as f:
            f.write(pdb_block)

class MoleculeSystem:

    def __init__(self,molecule_list:list[Molecule],molecule_num :list[int],density :float,name:str,tol=2.0,nloop=100):
        self.molecule_list = molecule_list
        self.molecule_num = molecule_num
        self.name = name
        self.__calc_length(density)
        self.__write_inp(tol=tol,nloop=nloop)

    def __calc_length (self,density):
        total_mol_wt = 0
        for mol,num in zip(self.molecule_list,self.molecule_num):
            total_mol_wt += mol.gram*num
        length_cm = (total_mol_wt / density)**(1/3)
        length_aa = length_cm * 1e8
        self.length_aa = length_aa

    def __write_inp(self,tol=2.0,nloop=100):
        packmol_input = f"tolerance {tol}\n"
        packmol_input += f"output {self.name}.pdb\n"
        packmol_input += f"pbc {self.length_aa} {self.length_aa} {self.length_aa}\n"
        packmol_input += "filetype pdb\n"
        for mol,num in zip(self.molecule_list,self.molecule_num):
            packmol_input += f"structure {mol.name}.pdb\n"
            packmol_input += f"  nloop {nloop}\n"
            packmol_input += f"  number {num}\n"
            packmol_input += f"  inside cube 0. 0. 0. {self.length_aa}\n"
            packmol_input += "end structure\n"
        with open(f"{self.name}.inp","w") as f:
            f.write(packmol_input)
        subprocess.call(f"{packmol_path} < {self.name}.inp",shell=True)

class PCFFSetup:

    def __init__(self,pdb_name):

        cur = os.path.abspath(os.curdir)
        atom_typing = {
            "-t": f"{os.path.abspath(pdb_name)}",
            "-ff": "PCFF"
        }
        all2_lmp = {
            "-t": f"{os.path.splitext(os.path.abspath(pdb_name))[0]+"_typed.data"}",
            "-n": f"{os.path.splitext(os.path.abspath(pdb_name))[0]+"_typed.nta"}",
            "-frc": "frc_files/pcff.frc"
        }
        os.chdir("/home/ryo/projects/LUNAR")
        command = "python atom_typing.py "
        command += " ".join([f"{key} {value}" for key, value in atom_typing.items()])
        subprocess.call(command,shell=True)
        command = "python all2lmp.py "
        command += " ".join([f"{key} {value}" for key, value in all2_lmp.items()])
        subprocess.call(command,shell=True)
        os.chdir(cur)

class PCFFInput:

    def __init__(self,datafile,cut=12):
        commands = []
        commands.append("units real")
        commands.append("atom_style full")
        commands.append("boundary p p p")
        commands.append(f"pair_style lj/class2/coul/long {cut}")
        commands.append("bond_style      class2")
        commands.append("angle_style     class2")
        commands.append("dihedral_style  class2")
        commands.append("improper_style  class2")
        commands.append("kspace_style    pppm 1.0e-4")
        commands.append("neighbor        2.0 bin")
        commands.append("neigh_modify    every 1 delay 0 check yes")
        commands.append(f"read_data {datafile}")
        commands.append("thermo_style custom step cpu cpuremain density etotal")
        commands.append("thermo          1000")
        commands.append("fix del_linear all momentum 1 linear 1 1 1 angular")
        self.commands = commands
        self.datafile = datafile
        self.count = 0

    def minimize(self,maxiter=10000):
        self.commands.append(f"minimize 1e-4 1e-6 {maxiter} 10000")
        self.commands.append(f"write_data {os.path.splitext(self.datafile)[0]}.0.data")

    def nvt(self,temp=300,tstop=300,tdamp=100,run_ns=1,xtc=True):
        self.count += 1
        self.commands.append(f"timestep 1.0")
        self.commands.append(f"fix {self.count} all nvt temp {temp} {tstop} {tdamp}")
        if xtc:
            self.commands.append(f"dump xtc{self.count} all atom 1000 {os.path.splitext(self.datafile)[0]}.{self.count}.lmptrj")
        self.commands.append(f"run {int(run_ns*1e6)}")
        self.commands.append(f"write_data {os.path.splitext(self.datafile)[0]}.{self.count}.data")
        self.commands.append(f"unfix {self.count}")

    def npt(self,temp=300,tstop=300,tdamp=100,p=1,pdamp=1000,run_ns=1,xtc=True):
        self.count += 1
        self.commands.append(f"timestep 1.0")
        self.commands.append(f"fix {self.count} all npt temp {temp} {tstop} {tdamp} iso {p} {p} {pdamp}")
        if xtc:
            self.commands.append(f"dump xtc{self.count} all atom 1000 {os.path.splitext(self.datafile)[0]}.{self.count}.lmptrj")
        self.commands.append(f"run {int(run_ns*1e6)}")
        self.commands.append(f"write_data {os.path.splitext(self.datafile)[0]}.{self.count}.data")
        self.commands.append(f"unfix {self.count}")

    def msd(self,query):
        self.commands.append(f"group group {query}")
        self.commands.append(f"compute msd group msd")
        self.commands.append(f"fix msdAve all ave/time 100 10 1000 c_msd file {os.path.splitext(self.datafile)[0]}.{self.count}.msd.txt")
        self.commands.append(f"unfix msdAve")

    def run(self,core=48,node=1,walltime="72:00:00"):
        print(os.path.abspath(os.path.splitext(self.datafile)[0]+".in"))
        with open(os.path.abspath(os.path.splitext(self.datafile)[0]+".in"),"w") as f:
            f.write("\n".join(self.commands))
        with open(os.path.abspath(os.path.splitext(self.datafile)[0]+".sh"),"w") as f:
            f.write("#!/bin/bash\n")
            f.write("#PBS -N lammps_job           # ジョブ名の指定\n")
            f.write(f"#PBS -l nodes={node}:ppn={core}       # 2ノード、各ノード16プロセッサ（例：合計32プロセス）\n")
            f.write(f"#PBS -l walltime={walltime}    # 最大実行時間（例：12時間）\n")
            f.write("#PBS -j oe                 # 標準出力とエラーを同じファイルに結合\n")
            f.write("#PBS -M your.email@example.com   # 通知メールの送信先\n")
            f.write("#PBS -m abe                # ジョブの開始、終了、異常終了時にメール送信\n")
            f.write("\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("module load lammps/15Dec2023\n")
            f.write(f"mpirun -np 32 lmp_mpi -in {os.path.splitext(self.datafile)[0]}.in\n")
        # ジョブディレクトリに移動（PBSで実行する場合は必須）
        # subprocess.call(f"lmp_gpu -sf hybrid gpu omp -pk gpu 1 omp 16 -in {os.path.splitext(self.datafile)[0]}.in",shell=True)

if __name__ == "__main__":
    os.chdir("/home/ryo/projects/LUNAR/TEST/")
    soln = Molecule(Chem.MolFromSmiles("[Si]1(C)(C)O"+"[Si](C)(C)O"*6+"[Si](C)(C)O1"),name="cODMS8")
    soln.mol_obj = Chem.AddHs(soln.mol_obj)
    solv = Molecule(Chem.MolFromSmiles("C1CCCCC1[Br]"),name="BrCH")
    solv.mol_obj = Chem.AddHs(solv.mol_obj)
    soln.load_pdb()
    solv.load_pdb()
    N_soln = 20
    wf = 5e-2
    N_solv = N_soln*soln.gram * (1-wf) / wf / solv.gram
    N_solv = int(N_solv) 
    system = MoleculeSystem([soln,solv],molecule_num=[N_soln,N_solv],density=0.7,name=f"{soln.name}_{solv.name}_pcff")
    PCFFSetup(f"../TEST/{soln.name}_{solv.name}_pcff.pdb")
    pcff_input = PCFFInput(f"/home/ryo/projects/LUNAR/TEST/{soln.name}_{solv.name}_pcff_typed_IFF.data")
    pcff_input.minimize(maxiter=50000)
    pcff_input.nvt(temp=300,tstop=300,run_ns=1)
    pcff_input.npt(run_ns=10,xtc=False)
    pcff_input.npt(run_ns=10,xtc=True)
    pcff_input.run()
