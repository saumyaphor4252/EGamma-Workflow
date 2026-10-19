
import subprocess

runs = [
    403751, 403774, 403775, 403811,
    403818, 403867, 403894, 403937
]

dataset_template = (
    "/EGamma{}/Run2026D-ZElectron-PromptReco-v1/RAW-RECO"
)

output_file = "das_files.txt"
total_files = 0

with open(output_file, "w") as out:
    for x in range(6): # EGamma{0,1,2,3,4,5} Datasets
        dataset = dataset_template.format(x)
        print(f"\nSearching dataset: {dataset}")

        for run in runs:
            query = f'file dataset={dataset} run={run}'

            result = subprocess.run(
                ["dasgoclient", "-query", query],
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                print(f"Error querying run {run}: {result.stderr.strip()}")
                continue

            files = result.stdout.strip().splitlines()

            for filename in files:
                filename = filename.strip()

                if filename:
                    out.write(f"file:/eos/cms/{filename}\n")
                    total_files += 1

            print(f"Run {run}: found {len(files)} files")

print(f"\nTotal files found: {total_files}")
print(f"File list saved to: {output_file}")
