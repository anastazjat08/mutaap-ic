import requests
import os
import time

def predictESM(version: str, sequence: str, out_dir: str) -> str:
    '''Predicts the structure of a protein sequence using ESMFold and saves the predicted structure as a PDB file.
    Parameters
    ----------
    version : str
        A string identifier for the version of the sequence being predicted ("orig", "mut").
    sequence : str
        The amino acid sequence of the protein to predict.
    out_dir : str
        The directory where the predicted structure will be saved.
    '''
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    print(f"Starting prediction for sequence: {sequence}")

    max_retries = 5
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.post(url, data=sequence, timeout=60)

            if response.status_code == 200:
                print("Successfully predicted structure. Saving to file...")
                # ensure output directory exists
                os.makedirs(out_dir, exist_ok=True)
                file_path = os.path.join(out_dir, f"{version}_esmfold_v1.pdb")
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(response.text)
                return file_path

            else:
                print(f"[Attempt {attempt}] API error {response.status_code}: {response.text}")

                # 504 = timeout on server side, so retry
                if response.status_code == 504:
                    time.sleep(5)
                    continue
                else:
                    break

        except requests.exceptions.Timeout:
            print(f"[Attempt {attempt}] Request timed out. Retrying...")
            time.sleep(5)

        except requests.exceptions.RequestException as e:
            print(f"[Attempt {attempt}] Network error: {e}")
            time.sleep(5)

    print(f"[ERROR] Failed to predict structure for {version} after {max_retries} attempts.")
    return None


 # TEST
# test_seq = "MKTAYIAKQRQISFVKSHFSRQDILDLIYQYARVVYQ"
# output_file = "/home/nastka/ADP/mutaap_test/predicted_structure.pdb"

# predicted_pdb = predictESM('orig', test_seq, "/home/nastka/ADP/mutaap_test")