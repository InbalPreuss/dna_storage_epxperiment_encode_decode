# DnaStorage

## Usage

Create a venv (recommended)

```console
git clone https://github.com/InbalPreuss/DnaStorage.git
pip install -r requirements.txt
```

inside DnaStorage directory, next to dna_storage directory create data/testing/input_text.dna with some text inside.

All the outputs will be written next to it.

then just run:
```console
python main.py
```

## Test
```console
python -m pytest
```
Create graphs:
In test_dna.py run test_number_of_oligos_per_barcode()

## Profiling
in test_dna.py make sure that the function in
```console
if __name__ == '__main__':
```
is:
```console
code_profiling()
```
Then, run:

linux:
```console
python -m cProfile -s time test_dna.py > temp_file && head -n 10 temp_file
```
Windows: 
```console
python -m cProfile -s time test_dna.py > profiling_data_no_synthsis_1_KB.txt
```


## TODO:
4. Explain in config that the barcode len is for ACGT, and the payload len is for Z1,Z2,... 3-Mers, 3 nucleotides
1. Add Multi threading
6. Check the time in the text handling - If it takes too long, we should remove this stage and generate random binary
4. rename results_file to output_file
5. add typings
6. add docstring
7. think of adding check after payload reed solomon


## Inbal:

1. Run KB Data and measure the time
2. Run MB Data and measure the time

## Done:

2. Make RS for all oligos
1. Make sure when adding and removing a letter, of reading a 3 letter sequence that does not exist, we would remove that specific oligo
2. The error graph should be logarithmic graph
3. Remove the exactly 5 bins to distribute the oligos in the synthesis phase - in real life it could be that not all oligos are chosen to In big numbers it probably wouldn't happen). We can cahnge the number of generated oligos to be 100 and not 20.
4. Make graphs for different number of oligos copies: 20, 100, 1000, 10000. Each graph should use all the M samples
5. Make RS for each oligo
6. After synthesising the DNA we should 
    1) Shuffle the data and then use the sort algorithm to sort the oligo according the barcode
7. Code that generates different data sizes
8. Make tests for RS
7. Remove x00 from text result.
2. After synthesising the DNA we should 
    1) Done! Shuffle the data and then use the sort algorithm to sort the oligo according the barcode
    2) Take M(= density of the sequencing) samples each time for the reading   
3. Make graphs for different M: 20, 50, 100. Each graph should use 10000 copies for one oligo.
3. Make RS for 16c3 
4. Make RS for 16c7
5. Put everything in packages -> test package, algorithm package...
