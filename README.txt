#Download Error_prediction.tar file
       tar -xvf Error_prediction.tar
       Rscript error_prediction.r  test.fasta   output.txt  eukaryote_genome.txt #to predict errors in eukaryote dataset
       Rscript error_prediction.r  test.fasta   output.txt  prokaryote_genome.txt #to predict errors in prokaryote dataset
	   
Note : Place "test.fasta" in the "Error_prediction" directory.
       test.fasta - name of the input dataset  provided by the user.
       output.txt - name of the output file.

#Files generated : 
- feature_matrix.txt : Include the features set[read_length,distance,ratio,nucleotide probabilities], generated from the input dataset provided by the User. 
- attribute_importance_matrix.txt : attribute imporatnce of each feature at lamda.min cutoff. 
- output.txt : Error predictions in the reads of Input dataset, provided by the user.		
		
#Test Run 
cd testrun
Rscript error_prediction_test.r  test.fasta   output.txt  prokaryote_genome.txt
		
		