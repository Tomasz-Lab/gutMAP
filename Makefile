# The gutMAP workflow on Helios: download, batch runs and promote.
# Run make in the launch dir, the dir that holds src/. Put the settings in config.mk
# or give them on the command line. `make help` lists the targets.
SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help

CONFIG ?= config.mk
-include $(CONFIG)

JOB_LOG_DIR ?= logs
LOG_DIR := $(abspath $(JOB_LOG_DIR))
SBATCH = sbatch --parsable -o $(LOG_DIR)/%x_%j.out $(if $(ACCOUNT),-A $(ACCOUNT)) $(if $(PARTITION),-p $(PARTITION))
DOWNLOAD_ARGS = --ids $(abspath $(IDS)) --input-dir $(abspath $(INPUT_DIR))
BATCH_ARGS = --staging-dir $(abspath $(STAGING_DIR)) $(DOWNLOAD_ARGS)
ME := $(shell id -un)
export KF_ENV

need = $(foreach v,$(1),$(if $($(v)),,$(error $(v) is not set. Set it in $(CONFIG) or on the command line)))
no_job = jobs=$$(squeue -h -u $(ME) -n $(1)); if [ -n "$$jobs" ]; then echo "A $(1) job is in the queue. $(2)" >&2; exit 1; fi

.PHONY: help runinfo download-plan download validate download-status plan run stop status sample promote queue collect-hjh

help: ## list the targets
	@grep -hE '^[a-z-]+:.*## ' $(firstword $(MAKEFILE_LIST)) | sed -E 's/:.*## /\t/' | expand -t 18

runinfo: ## add the NCBI runinfo of new samples to <INPUT_DIR>.download/runinfo.csv
	@$(call need,IDS INPUT_DIR)
	python3 src/download/download.py runinfo $(DOWNLOAD_ARGS)

download-plan: ## show which samples the download job will get
	@$(call need,IDS INPUT_DIR)
	python3 src/download/download.py plan $(DOWNLOAD_ARGS)

download: ## submit the download job and the validate job after it
	@$(call need,IDS INPUT_DIR KF_ENV)
	@$(call no_job,gutmap-batches,Downloads and batches must not run at the same time.)
	mkdir -p $(LOG_DIR)
	cd src/download && J=$$($(SBATCH) download.sbatch $(DOWNLOAD_ARGS) $(if $(METHODS),--methods $(METHODS))) \
	    && echo "download job $$J" \
	    && V=$$($(SBATCH) --dependency=singleton,afterany:$$J validate.sbatch --input-dir $(abspath $(INPUT_DIR))) \
	    && echo "validate job $$V"

validate: ## submit only the validate job
	@$(call need,INPUT_DIR KF_ENV)
	mkdir -p $(LOG_DIR)
	cd src/download && V=$$($(SBATCH) validate.sbatch --input-dir $(abspath $(INPUT_DIR))) && echo "validate job $$V"

download-status: ## summary of the download attempts and content checks
	@$(call need,INPUT_DIR)
	python3 src/download/download.py status --input-dir $(abspath $(INPUT_DIR))

plan: ## show the next batch
	@$(call need,STAGING_DIR IDS INPUT_DIR DISK_BUDGET)
	python3 src/batches/batches.py plan $(BATCH_ARGS) --disk-budget $(DISK_BUDGET)

run: ## submit the batches job (N_BATCHES=n to stop after n batches)
	@$(call need,STAGING_DIR IDS INPUT_DIR DISK_BUDGET)
	@$(call no_job,gutmap-download,Downloads and batches must not run at the same time.)
	@$(call no_job,gutmap-validate,Downloads and batches must not run at the same time.)
	mkdir -p $(LOG_DIR)
	J=$$($(SBATCH) src/batches/batches.sbatch $(BATCH_ARGS) --disk-budget $(DISK_BUDGET) $(if $(N_BATCHES),--n-batches $(N_BATCHES))) \
	    && echo "batches job $$J"

stop: ## stop the batches job; the next run finishes its batch
	scancel --name=gutmap-batches --user=$(ME)

status: ## sample states and the result of each batch
	@$(call need,STAGING_DIR IDS INPUT_DIR DISK_BUDGET)
	python3 src/batches/batches.py status $(BATCH_ARGS) --disk-budget $(DISK_BUDGET)

sample: ## the row of one sample (SAMPLE=SRR123)
	@$(call need,STAGING_DIR SAMPLE)
	python3 src/batches/batches.py sample --staging-dir $(abspath $(STAGING_DIR)) $(SAMPLE)

promote: ## move the accepted samples (dryrun=1 to only show the counts)
	@$(call need,STAGING_DIR VALID_DIR VALID_LOGS_DIR)
	python3 src/batches/batches.py promote --staging-dir $(abspath $(STAGING_DIR)) --valid-dir $(abspath $(VALID_DIR)) \
	    --logs-dir $(abspath $(VALID_LOGS_DIR)) $(if $(dryrun),--dry-run)

queue: ## the Slurm jobs of this user
	squeue -u $(ME) -o "%.10i %.32j %.8T %.10M %R"

collect-hjh: ## add the Slurm job history with billing to HJH_CSV
	@$(call need,HJH_CSV)
	python3 utils/collect_hjh.py --out $(abspath $(HJH_CSV)) $(if $(HJH_START),--start $(HJH_START))
