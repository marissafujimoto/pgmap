import csv


def read_barcodes(barcodes_path: str) -> dict[str, str]:
    barcode_mapping = {}

    with open(barcodes_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        
        for barcode, sample_id in tsv_reader:
            barcode_mapping[barcode] = sample_id
    
    return barcode_mapping