Human Gut Microbiome project

## Run the workflow on Helios

The `Makefile` runs the workflow from the launch dir, the dir that holds `src/`.

1. Set the paths, the disk budget, the conda env of kingfisher and the Slurm account in `config.mk`.
2. Get the runinfo and download the samples:

   ```
   make runinfo
   make download-plan
   make download
   make download-status
   ```

3. When the download and validate jobs end, run the batches:

   ```
   make plan
   make run
   make status
   ```

4. Promote the accepted samples. Without `YES=1`, `promote` only shows the counts:

   ```
   make promote
   make promote YES=1
   ```

- `make stop` stops the batches job after it finishes the batch that runs. `make run` continues the work.
- `make sample SAMPLE=SRR123` shows the state of one sample and of each of its steps.
- `make queue` shows your Slurm jobs. `make help` lists all targets.
- Downloads and batches must not run at the same time. `make download` and `make run` stop if a job of the other kind is in the queue.
- A setting on the command line replaces the one in `config.mk`, for example `make run N_BATCHES=1`.

For the details, see `src/download/README.md` and `src/batches/README.md`.
