import subprocess
import sys

# Specific function test
sys.stderr.write("Running specific function tests\n")
subprocess.run("python ../tape_specific-runner.py --spec-amino 50 --distance 50 --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 50 --distance 50 --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 50 --distance 100 --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("N-terminus specific function test\n")
subprocess.run("python ../tape_specific-runner.py --spec-amino 1 --distance 50 --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 1 --distance 50 --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 1 --distance 50 --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("C-terminus specific function test\n")
subprocess.run("python ../tape_specific-runner.py --spec-amino 3739 --distance 50 --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 3739 --distance 50 --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 3739 --distance 50 --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

# General function test
sys.stderr.write("Running general function tests\n")
subprocess.run("python ../tape_general-runner.py --aa K --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa K --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa K --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("N-terminus general function test\n")
subprocess.run("python ../tape_general-runner.py --aa M --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa M --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa M --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("C-terminus general function test\n")
subprocess.run("python ../tape_general-runner.py --aa - --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa - --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa - --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("Testing multi-aa targeting\n")
subprocess.run("python ../tape_general-runner.py --aa KL --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa KL --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa KL --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("N-terminus multi-aa general function test\n")
subprocess.run("python ../tape_general-runner.py --aa MV --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa MV --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa MV --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("C-terminus multi-aa general function test\n")
subprocess.run("python ../tape_general-runner.py --aa A- --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa A- --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa A- --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
