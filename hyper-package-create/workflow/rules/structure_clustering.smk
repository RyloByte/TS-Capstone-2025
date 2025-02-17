import itertools
import os

PDB_DIR = "data/pdb_structures"
TMSCORE_EXEC = "./TMscore"

pdb_files = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
pairs = [(os.path.join(PDB_DIR, f1), os.path.join(PDB_DIR, f2)) for f1, f2 in itertools.combinations(pdb_files, 2)]

rule pairwise_tm_score:
    input:
        pdb1=lambda wildcards: wildcards.pdb1,
        pdb2=lambda wildcards: wildcards.pdb2
    output:
        "data/{pdb1.stem}_vs_{pdb2.stem}.txt"
    shell:
        """
        {TMSCORE_EXEC} {input.pdb1} {input.pdb2} | grep 'TM-score    =' | awk '{print $3}' > {output}
        """

rule aggregate_tm_scores:
    input:
        expand("results/{{pdb1.stem}}_vs_{{pdb2.stem}}.txt", pdb1=[p.split("/")[-1] for p in pdb_files], pdb2=[p.split("/")[-1] for p in pdb_files])
    output:
        "data/tm_results.csv"
    run:
        with open("data/tm_results.csv", "w") as out_file:
            out_file.write("Structure1,Structure2,TM-score\n")
            for result_file in input:
                score_line = f.readline().strip().split("=")
                if len(score_line) > 1:
                    score = score_line[1]
                    pdb1, pdb2 = os.path.basename(result_file).replace(".txt", "").split("_vs_")
                    outfile.write(f"{pdb1},{pdb2},{score}\n")

# Cluster based on tm scores, split fastas based on that
