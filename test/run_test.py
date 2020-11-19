import subprocess
import sys

# Specific function test
sys.stderr.write("Running specific function tests\n")
subprocess.run("python ../tape_specific-runner.py --spec-amino 50 --distance 100 --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 50 --distance 100 --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 50 --distance 100 --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("N-terminus specific function test\n")
subprocess.run("python ../tape_specific-runner.py --spec-amino 1 --distance 100 --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 1 --distance 100 --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --spec-amino 1 --distance 100 --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

sys.stderr.write("C-terminus specific function test\n")
subprocess.run("python ../tape_specific-runner.py --aa-pos 50 --distance 100 --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --aa-pos 50 --distance 100 --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_specific-runner.py --aa-pos 50 --distance 100 --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)

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
subprocess.run("python ../tape_general-runner.py --aa * --motif NGG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa * --motif YG --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
subprocess.run("python ../tape_general-runner.py --aa * --motif TTTN --loci test_loci.txt --cds test_cds.txt --genome test_genome.txt --output test", shell=True, check=True)
