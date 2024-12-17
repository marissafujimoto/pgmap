import csv


def read_barcodes(barcodes_path: str) -> dict[str, str]:
    """
    Read barcodes from a barcode tsv file. The first column of the tsv should be the barcode sequences and the
    second should be the sample id corresponding to that barcode.

    Args:
        barcodes_path (str): The path to a barcodes file (a tsv).

    Returns:
        barcodes_mapping (dict[str, str]): A dictionary mapping from each barcode to the corresponding sample id.
    """
    barcode_mapping = {}

    with open(barcodes_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        
        for barcode, sample_id in tsv_reader:
            barcode_mapping[barcode] = sample_id
    
    return barcode_mapping