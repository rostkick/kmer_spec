# kmer_spec
### How to install it
You can clone whole repo by standart command:
```
git clone https://github.com/rostkick/kmer_spec.git
```
Also if you need to save it **without commit-log** you schould use:
```
git clone —depth=1 https://github.com/rostkick/kmer_spec.git
```
If you want to save **only one** file (for example without saving small_test.fastq) you should:
  1) copy link of favorite file
  2) go to https://minhaskamal.github.io/DownGit/#/home, paste it, and create Download Link.
### Keys info
  
  ```-i, --input <str>``` Paste your path to input .fastq file here.  
  ```-k, --kmer <int>``` Set k-mer length (For example 23)  
  ```-q, --quality <int>``` Set reads quality threshold (For example 30)  
  ```-w, --write <bool>``` You want to save the plot (.png)?   
  ```-xmin, --borderxmin <int>``` Set graphic borders (xmin)  
  ```-xmax, --borderxmax <int>``` Set graphic borders (xmax)  
  ```-ymin, --borderymin <int>``` Set graphic borders (ymin)  
  ```-ymax, --borderymax <int>``` Set graphic borders (ymax)  
  ```-ax, --axis <str>``` If you want use auto borders -- set 'auto', else -- don't use it.  
  ```-gs, --get_spectr <bool>``` Write .tsv file with spectrum for visualize somewhere else'
### How to use it
```
python3 kmer_spectr.py -i small_test.fastq -k 15 -q 20 -w True
```
![alt text](https://github.com/rostkick/kmer_spec/blob/master/small_test.fq_k15-q20.png)

More life example
```
python3 kmer_spectr.py -i test_kmer.fastq -k 15 -q 20 -w True
```
![alt text](https://github.com/rostkick/kmer_spec/blob/master/test_kmer.fastq_k15-q20.png)

If you want just get k-mer spectrum withou visualization use it:
```
python3 kmer_spect -i test_kmer.fastq -k 31 -q 20 --get_spectr True -w False
```
### Noise cut-off and genome size calculating.
Finding minimum in small distance between 2 first picks. The cut-off is calculated automatically.
Then the approximated size of the genome will be calculating with depending of the cut-off.
