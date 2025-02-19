This program is a C++ implementation of an N-body simulation. 

But before you begin, make sure you have access to:
- A Linux-based system with SLURM workload management
- gcc compiler (this is critical to execute the program)
- Access to hpc lab computer Centaurus
- UNC Charlotte VPN if you are outside of campus or not using "eduroam" network
- Python 3 with pandas and matplotlib libraries for plot.py


Steps to compile and experiment:
1. Connecting hpc lab computer by "ssh hpc-student.charlotte.edu -l your-username"
2. Authenticate with Duo
3. Type "g++ A2.cpp -o A2" A2 is the name of the executable file and A2.cpp is the source code. g++ allows us to make executable file that is binary so CPU can process the high level program.
4. Schedule the job by "sbatch A2.sh"
5. Outcome should be something like "Submitted batch job [????]", pay attention to the number.
6. Wait a bit for command to finish running and record the time it takes.
7. There are 3 jobs scheduled in .sh file. To look at the time it took for each process, open record.txt by typing "cat record.txt"
8. To create a chart , type "python3 plot.py solar.tsv solar.pdf 10000" This will create pdf file of plot data so you can see  visually.


Keys:
- If  you would like to run the progarm with your sets of data, you can execute it by "./A1 <number of elements(particles)> <Time  step> < Number of time steps> <Output state interval>"

Records from my  end:
- 100 1 10000 100, Run time was: 2155 milliseconds.
- 2 200 5000000 100,  Run time was: 1356 milliseconds.
- 1000 1 10000 100, Run time was: 206766 milliseconds.


