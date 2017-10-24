% GB comments
1.	100
2a. 70 Close but the totatlen is calculated using only a fraction of the total sequence. Should have used seq1 or seq2 to determine the entire length. 
2b. 70 same issue as 2a
2c. 70 same issue as 2a. 
3a 100 
3b. 100
3c. 100  	
Overall: 87


%HW3

%% Problem 1 - Smith-Waterman alignment
%Walter Frank Lenoir
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

%Please see "wflenoir-hw3-1.JPG" Optical alignment has a result of 10, or 
%11 if the end gap is not included. 

%% Problem 2 - using the NCBI databases and sequence alignments
%Walter Frank Lenoir

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

ERK1 = getgenbank('NM_002746');
ERK2 = getgenbank('NM_002745');
index1 = ERK1.CDS.indices;
index2 = ERK2.CDS.indices;
seq1 = ERK1.Sequence(index1(1):index1(2));
seq2 = ERK2.Sequence(index2(1):index2(2));

[score,align,start] = swalign(seq1,seq2,'Alphabet','nt');
matches = count(align(2,:),'|');
totlen = length(align);

%The fraction of basepairs that in ERK1 that can align to ERK2 is 811/1073
%(matches/totlen).

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

aaseq1 = ERK1.CDS.translation;
aaseq2 = ERK2.CDS.translation;

[score,align,start] = swalign(aaseq1,aaseq2);
matches = count(align(2,:),'|');
totlen = length(align);

showalignment(align);

%There fraction of amino acids that align is 305/346. 96% of the alignment had positive matches (functional
%matching of amino acids)

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

musERK1 = getgenbank('NM_011952'); %mus musculus erk1 mrna
musERK2 = getgenbank('NM_011949'); %mus musculus erk2 mrna
index1 = musERK1.CDS.indices;
index2 = musERK2.CDS.indices;
musseq1 = musERK1.Sequence(index1(1):index1(2));
musseq2 = musERK2.Sequence(index2(1):index2(2));
musaaseq1 = musERK1.CDS.translation;
musaaseq2 = musERK2.CDS.translation;

[score,align,start] = swalign(seq1,musseq1,'Alphabet','nt');
matches = count(align(2,:),'|');
totlen = length(align);

%For ERK1, the human and mouse sequences had 1026 matches/1137 total length
%in the alignment

[score,align,start] = swalign(seq2,musseq2,'Alphabet','nt');
matches = count(align(2,:),'|');
totlen = length(align);

%For ERK2, the human and mouse sequences had 996 matches /1075 total length
%in the alignment.

[score,align,start] = swalign(aaseq1,musaaseq1);
matches = count(align(2,:),'|');
totlen = length(align);

showalignment(align);


%For ERK1, the human and mouse amino acids had 367 matches /378 total length
%in the alignment. 98% of the alignment had positive matches (functional
%matching of amino acids).


[score,align,start] = swalign(aaseq2,musaaseq2);
matches = count(align(2,:),'|');
totlen = length(align);

showalignment(align);

%For ERK2, the human and mouse amino acids had 355 matches /357 total length
%in the alignment. 100% of the alignment had positive matches (functional
%matching of amino acids).

%Overall the amino acid sequences were very similar. ERK2 appears to be
%more similar in terms of translational sequence (99% vs 97% for ERK1), and 
% in terms of nucleotide sequence (92.5% vs 90% for ERK1). 

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

result = blasthits('NM_011949',50); % can only take 50, max hits parameter in blast does not appear to work. 

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

[humanhit,organismhit] = orgcompare('NM_011949'); %warnings will come up but no errors

%humanhit = NM_002745
%organismhit = XM_021184651

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

[humanhit,organismhit] = orgcompare('NM_001126118'); %TP53

%Humanhit = NM_001126112.2
%Organismhit = XM_016931470.1
%The human hit is another version of TP53 in ncbi, while the closest animal
%hit is TP53 in chimpanzees. TP53 likely has multiple copies in NCBI. These
%results make sense. 

[humanhit,organismhit] = orgcompare('EU787372'); 

%This is a gene that comes from the Emys orbicularis (European pond
%turtle). There was no human hit, however the closest organism hit was
%'XM_007063781.1', a gene from the green sea turtle. There was likely no
%human hit because blast returns only 50 hits (max hit return does not appear
%to be a variable type anymore), and there were at least 50 other genes more 
%closely related from other species than humans. There are at least 50 
%other gene hits in other species between the examined gene and humans. 

