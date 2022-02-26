import itertools
from textwrap import wrap
from typing import Union, Dict, List, Tuple
from pathlib import Path

import pandas as pd

from dna_storage.rs_adapter import RSBarcodeAdapter, RSPayloadAdapter, RSWideAdapter
from dna_storage import utils


#################################################################
# @ Class: Encoder
# @ Description: Retrieve the oligo to the oligo that was written
#                in originally
#################################################################

class Encoder:
    def __init__(self, barcode_len: int,
                 barcode_rs_len: int,
                 payload_len: int,
                 payload_rs_len: int,
                 binary_file_name: str,
                 shrink_dict: Dict,
                 k_mer: int,
                 k_mer_representative_to_z: Dict,
                 binary_to_z: Dict,
                 subset_size: int,
                 bits_per_z: int,
                 oligos_per_block_len: int,
                 oligos_per_block_rs_len: int,
                 barcode_coder: RSBarcodeAdapter,
                 payload_coder: RSPayloadAdapter,
                 wide_coder: RSWideAdapter,
                 results_file: Union[Path, str],
                 results_file_without_rs_wide: Union[Path, str],
                 barcode_dict: Dict,
                 z_to_k_mer_representative: Dict
                 ):
        self.file_name = binary_file_name
        self.barcode_len = barcode_len
        self.barcode_rs_len = barcode_rs_len
        self.payload_len = payload_len
        self.payload_rs_len = payload_rs_len
        self.shrink_dict = shrink_dict
        self.k_mer = k_mer
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.binary_to_z_dict = binary_to_z
        self.subset_size = subset_size
        self.bits_per_z = bits_per_z
        self.oligos_per_block_len = oligos_per_block_len
        self.oligos_per_block_rs_len = oligos_per_block_rs_len
        self.results_file = results_file
        self.results_file_without_rs_wide = results_file_without_rs_wide
        open(self.results_file, 'w').close()
        self.barcode_generator = utils.dna_sequence_generator(sequence_len=self.barcode_len)
        self.barcode_generator_pre_determend_barcode = utils.dna_sequence_generator_pre_determend_barcode(barcode_dict=barcode_dict)
        self.barcode_coder = barcode_coder
        self.payload_coder = payload_coder
        self.wide_coder = wide_coder
        self.z_to_k_mer_representative = z_to_k_mer_representative

    def run(self):
        number_of_blocks = 0
        with open(self.file_name, 'r', encoding='utf-8') as file:
            z_list_accumulation_per_block = []
            for line in file:
                line = line.strip('\n')
                z_list = []
                for binary_to_transform in wrap(line, self.bits_per_z):
                    z = self.binary_to_z(binary=binary_to_transform)
                    z_list.append(z)
                z_list_accumulation_per_block.append(z_list)

            z_list_accumulation_with_rs = z_list_accumulation_per_block
            amount_oligos_per_block_len_to_write = 0
            for z_list in z_list_accumulation_with_rs:
                oligo = self.z_to_oligo(z_list)
                self.save_oligo(results_file=self.results_file, oligo=oligo)
                if amount_oligos_per_block_len_to_write > 0:
                    self.save_oligo(results_file=self.results_file_without_rs_wide, oligo=oligo)
                    amount_oligos_per_block_len_to_write = amount_oligos_per_block_len_to_write - 1

        self.create_source_target_well_for_robot()
        return number_of_blocks

    def binary_to_z(self, binary: str) -> str:
        binary_tuple = tuple([int(b) for b in binary])
        return self.binary_to_z_dict[binary_tuple]

    def z_to_oligo(self, z_list: List[str]) -> str:
        oligo = z_list
        barcode = next(self.barcode_generator)
        barcode = "".join(barcode)
        oligo.insert(0, barcode)
        return ",".join(oligo)

    def wide_block_rs(self, z_list_accumulation_per_block: List[List[str]]) -> List[List[str]]:
        rs_append = [[] for _ in range(int(self.oligos_per_block_len + self.oligos_per_block_rs_len))]
        for col in range(len(z_list_accumulation_per_block[0])):
            z_list = [elem[col] for elem in z_list_accumulation_per_block]
            col_with_rs = self.add_payload_rs_symbols_for_error_correction(payload=z_list, payload_or_wide='wide')
            for idx, z in enumerate(col_with_rs):
                rs_append[idx].append(z)
        return rs_append

    def add_payload_rs_symbols_for_error_correction(self, payload: Union[str, List[str]],
                                                    payload_or_wide: str = 'payload') -> List[str]:
        if isinstance(payload, str):
            payload = [c for c in payload]
        if payload_or_wide == 'payload':
            payload_encoded = self.payload_coder.encode(payload)
        else:
            payload_encoded = self.wide_coder.encode(payload)

        return payload_encoded

    def add_barcode_rs_symbols_for_error_correction(self, barcode: Tuple[str]) -> List[str]:
        barcode = list(barcode)
        barcode_encoded = self.barcode_coder.encode(barcode=barcode)
        return barcode_encoded

    def save_oligo(self, results_file: Union[Path, str], oligo: str) -> None:
        with open(results_file, 'a+', encoding='utf-8') as f:
            f.write(oligo + '\n')

    def create_source_target_well_for_robot(self):
        data = pd.read_csv(self.results_file, header=None)
        df = pd.DataFrame()
        df_list = []
        for _, row in data.iterrows():
            d = {}
            # d[0] = row[0]
            idx=1
            for i in range(1,5):
                z = row[i]
                xs = self.z_to_k_mer_representative[z]
                for j,x in enumerate(xs,1):
                    d[idx] = x
                    idx+=1
            df_list.append(d)

        final_df_list = []
        A = list(range(1,22))
        B = ["A","C","E","G","I","K","M","O"]
        symbols = [f"{y}{x}"for x in A for y in B]
        col_index = [0] * 5 + [1] * 5 + [2] * 5 + [3] * 5
        index = 0
        import random

        # for row in range(len(df_list)):
        #     for acc in range(0, 4):
        #         xs = []
        #         for j in range(1, 6):
        #             col = j + acc
        #             xs.append(df_list[row][col])
        #         random.shuffle(xs)
        #         for idx, j in enumerate(range(1, 6)):
        #             col = j + acc
        #             df_list[row][col] = xs[idx]

        # for row in range(len(df_list)):
        #     xs = list(df_list[row].values())
        #     random.shuffle(xs)
        #     for col in range(1, 21):
        #         df_list[row][col] = xs[col-1]

        for row in range(len(df_list)):
            for acc in range(0, 4):
                for j in range(1, 6):
                    col = j + acc * 5
                    print(df_list[row][col])
                    print(f"X{int(df_list[row][col][1:]) + acc * 16}")
                    df_list[row][col] = f"X{int(df_list[row][col][1:]) + acc * 16}"
            xs = list(df_list[row].values())
            random.shuffle(xs)
            for col in range(1, 21):
                df_list[row][col] = xs[col - 1]

        for i, col in zip(col_index, df_list[0].keys()):
            index=index%167

            for row in range(len(df_list)):
                number = int(df_list[row][col][1:])
                d = {"Source_Well": number, "Target_Well": symbols[index]}
                print(index)
                final_df_list.append(d)
                index+=1

        final_df = pd.DataFrame(final_df_list)
        final_df.to_excel("C:\\Users\\Inbal\\PycharmProjects\\DnaStorage_experiment_BarIlan_encode_decode\\data\\final_df.xlsx")

        a=3


