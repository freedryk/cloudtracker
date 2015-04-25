[Cloud-tracker](https://github.com/lorenghoh/cloud-tracker "cloud-tracker")
==============

## Updated cloud-tracking algorithm for SAM ##

## Example ##
To run the cloud-tracking algorithm, do
	
	python -O -B -m kernprof -l -v run_cloudtracker.py > line_stats.txt

or, to run with a memory profiler instead,	

	python -O -B -m memory_profiler run_cloudtracker.py > memory_stats.txt
