#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
from utilities import create_new_project, check_project_status, launch_jobsub, check_git_branch
from typing import Optional

MEDULLA_REPO_URL = "https://github.com/Fatimakhalifaa/medulla_numu.git"


def main(
    project_dir: str,
    experiment: str,
    create_project: bool,
    launch_jobs: int = None,
    test_job: bool = False,
    tml: str = None,
    batch_size: int = None,
    systematic: str = None,
    branch: str = "develop",
    memory: int = 1800,
    disk: Optional[int] = None,
    lifetime: str = "1h",
    relaunch_missing: bool = False,
):
    project_dir = Path(project_dir)
    project_exists = (project_dir / "project.db").exists()

    if create_project:
        if tml is None:
            raise ValueError("TOML file must be provided when creating a new project.")
        if batch_size is None:
            raise ValueError("Batch size must be provided when creating a new project.")
        if project_exists:
            raise FileExistsError(f"Project database {project_dir / 'project.db'} already exists.")

        create_new_project(project_dir, tml, batch_size, systematic)
        print(f"[INFO] -- Created new project in {project_dir}")

    if project_exists:
        check_project_status(project_dir)

    if test_job:
        if not project_exists:
            raise FileNotFoundError(
                f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first."
            )
        launch_jobsub(
            project_dir,
            experiment,
            njobs=1,
            branch=branch,
            memory=memory,
            disk=disk,
            lifetime=lifetime,
            relaunch_missing=relaunch_missing,
        )

    if launch_jobs is not None:
        if not project_exists:
            raise FileNotFoundError(
                f"Project database {project_dir / 'project.db'} does not exist. Please create a new project first."
            )
        launch_jobsub(
            project_dir,
            experiment,
            njobs=launch_jobs,
            branch=branch,
            memory=memory,
            disk=disk,
            lifetime=lifetime,
            relaunch_missing=relaunch_missing,
        )


if __name__ == "__main__":
    p = ArgumentParser(description="Run medulla.")

    p.add_argument("--project-dir", "-p", type=str, required=True)
    p.add_argument("--experiment", "-e", type=str, default="sbnd")
    p.add_argument("--create-project", "-c", action="store_true")
    p.add_argument("--toml", "-t", type=str)
    p.add_argument("--batch-size", "-b", type=int)
    p.add_argument("--systematic", "-s", type=str, default=None)
    p.add_argument("--test-job", "-T", action="store_true")
    p.add_argument("--launch-jobs", "-l", type=int, nargs="?", const=-1)
    p.add_argument("--relaunch-missing", "-r", action="store_true")
    p.add_argument("--branch", "-B", type=str, default="develop")
    p.add_argument("--memory", "-m", type=int, default=1800)
    p.add_argument("--disk", "-d", type=int, default=None)
    p.add_argument("--lifetime", "-f", type=str, default="1h")

    args = p.parse_args()

    if args.experiment not in ["sbnd", "icarus"]:
        p.error("Experiment must be either 'sbnd' or 'icarus'.")

    if args.create_project and args.toml is None:
        p.error("--toml is required when --create-project is set.")
    if args.create_project and args.batch_size is None:
        p.error("--batch-size is required when --create-project is set.")

    if args.test_job and args.launch_jobs is not None:
        p.error("--test-job and --launch-jobs are mutually exclusive.")

    if args.branch != "develop":
        print(f"[INFO] -- Using branch '{args.branch}' from fork repository.")
        if not check_git_branch(args.branch, MEDULLA_REPO_URL):
            p.error(
                f"Branch '{args.branch}' does not exist in {MEDULLA_REPO_URL}"
            )

    main(
        project_dir=args.project_dir,
        experiment=args.experiment,
        create_project=args.create_project,
        launch_jobs=args.launch_jobs,
        test_job=args.test_job,
        tml=args.toml,
        batch_size=args.batch_size,
        systematic=args.systematic,
        branch=args.branch,
        memory=args.memory,
        disk=args.disk,
        lifetime=args.lifetime,
        relaunch_missing=args.relaunch_missing,
    )
