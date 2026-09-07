"""Tests for ``publish_docs_gh-pages.sh``.

The script publishes a built site to a Pages branch as a single parentless commit,
force-pushed. Everything it does is exercised here against a bare repository used as
the remote, so the tests need no network and no GitHub credentials.

"""

import pathlib
import subprocess

import pytest

SCRIPT = pathlib.Path(__file__).parent.parent / "publish_docs_gh-pages.sh"

# Files written into the site directory by the `site_repo` fixture. The leading
# underscore and the nested directory are both cases the published branch must keep.
SITE_FILES = {
    "index.html": "<html>index</html>",
    "process_plate1.html": "<html>plate 1</html>",
    "sub/nested.html": "<html>nested</html>",
    "_leading_underscore.html": "<html>underscore</html>",
}


def git(*args, cwd):
    """Run a `git` command, returning its stripped stdout."""
    return subprocess.run(
        ["git", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def run_script(cwd, env_extra=None, script=SCRIPT):
    """Run the publish script, returning the `CompletedProcess` without checking it."""
    env = {
        "PATH": "/usr/bin:/bin:/usr/local/bin",
        # Isolate from the invoking user's git configuration, so the test does not
        # depend on their `core.excludesFile`, `init.defaultBranch`, or identity.
        "HOME": str(cwd),
        "GIT_CONFIG_NOSYSTEM": "1",
    }
    env.update(env_extra or {})
    return subprocess.run(
        ["bash", str(script)],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,  # the tests assert on the return code themselves
    )


@pytest.fixture
def site_repo(tmp_path):
    """A repo with a built site in `results/docs` and a bare repo as its remote."""
    remote = tmp_path / "remote.git"
    proj = tmp_path / "proj"
    git("init", "-q", "--bare", str(remote), cwd=tmp_path)
    git("init", "-q", str(proj), cwd=tmp_path)
    git("config", "user.name", "Test User", cwd=proj)
    git("config", "user.email", "test@example.com", cwd=proj)
    git("remote", "add", "origin", str(remote), cwd=proj)

    # The lab default `.gitignore`, which ignores the whole of `results/` and every
    # `_`-prefixed and dot-prefixed file. The published site must contain those files
    # regardless, so this is what would catch the site being staged with ignore rules
    # in force.
    (proj / ".gitignore").write_text(".*\n!.git*\n_*\n*.out\nresults/**\n")
    git("add", ".gitignore", cwd=proj)
    git("commit", "-q", "-m", "init", cwd=proj)
    git("push", "-q", "origin", "HEAD:refs/heads/main", cwd=proj)
    # `git ls-remote --symref` only reports a default branch if the remote's HEAD
    # resolves, which for a repo created by `git init --bare` it does not until set.
    git("symbolic-ref", "HEAD", "refs/heads/main", cwd=remote)

    write_site(proj, SITE_FILES)
    return proj, remote


def write_site(proj, files):
    """Write `files` (a path -> contents mapping) into the project's site directory."""
    site = proj / "results" / "docs"
    for name, contents in files.items():
        path = site / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)


def published_files(remote, branch="gh-pages"):
    """The set of file paths on `branch` of the remote."""
    listing = git("ls-tree", "-r", "--name-only", branch, cwd=remote)
    return set(listing.splitlines())


def test_publishes_exact_snapshot(site_repo):
    """The branch holds exactly the site's files, plus `.nojekyll`."""
    proj, remote = site_repo
    result = run_script(proj)
    assert result.returncode == 0, result.stderr
    assert published_files(remote) == set(SITE_FILES) | {".nojekyll"}
    for name, contents in SITE_FILES.items():
        assert git("cat-file", "-p", f"gh-pages:{name}", cwd=remote) == contents


def test_snapshot_commit_has_no_history(site_repo):
    """Publishing twice leaves one parentless commit, not a growing history."""
    proj, remote = site_repo
    assert run_script(proj).returncode == 0
    assert run_script(proj).returncode == 0
    assert git("rev-list", "--count", "gh-pages", cwd=remote) == "1"
    assert git("log", "-1", "--format=%P", "gh-pages", cwd=remote) == ""


def test_republishing_drops_deleted_files(site_repo):
    """A file removed from the site directory disappears from the branch."""
    proj, remote = site_repo
    assert run_script(proj).returncode == 0
    (proj / "results" / "docs" / "process_plate1.html").unlink()
    write_site(proj, {"index.html": "<html>index v2</html>"})

    assert run_script(proj).returncode == 0
    assert "process_plate1.html" not in published_files(remote)
    assert (
        git("cat-file", "-p", "gh-pages:index.html", cwd=remote)
        == "<html>index v2</html>"
    )


def test_leaves_no_trace_in_the_calling_repo(site_repo):
    """No local branch, staging ref, worktree, or scratch survives the publish."""
    proj, _remote = site_repo
    tmpdir = proj.parent / "tmpdir"
    tmpdir.mkdir()

    assert run_script(proj, {"TMPDIR": str(tmpdir)}).returncode == 0

    assert git("branch", "--list", "gh-pages", cwd=proj) == ""
    assert git("for-each-ref", "refs/publish-pages", cwd=proj) == ""
    assert git("worktree", "list", "--porcelain", cwd=proj).count("worktree ") == 1
    assert list((proj / ".git").glob("publish-pages.*")) == []
    # The script must not use TMPDIR at all: a site of any size needs no scratch space.
    assert list(tmpdir.iterdir()) == []


def test_does_not_disturb_another_runs_refs(site_repo):
    """A publish leaves refs under `refs/publish-pages/` alone.

    The snapshot commit is pushed by object name, so nothing is written under that
    namespace and nothing is cleaned up from it. Were a shared ref reintroduced there,
    together with a cleanup that deletes it unconditionally, two overlapping publishes
    could delete each other's in-flight commit.

    """
    proj, _remote = site_repo
    other = git("rev-parse", "HEAD", cwd=proj)
    git("update-ref", "refs/publish-pages/other-run", other, cwd=proj)

    assert run_script(proj).returncode == 0

    assert git("rev-parse", "refs/publish-pages/other-run", cwd=proj) == other


def test_force_staging_defeats_git_dir_excludes(site_repo):
    """`$GIT_DIR/info/exclude` cannot drop files out of the published site.

    Redirecting the work tree to the site directory already means the repository's
    top-level `.gitignore` is not consulted, but the git directory's own exclude file
    still is, so staging has to override it.

    """
    proj, remote = site_repo
    (proj / ".git" / "info" / "exclude").write_text("*\n")

    assert run_script(proj).returncode == 0
    assert published_files(remote) == set(SITE_FILES) | {".nojekyll"}


def test_refuses_to_publish_a_snapshot_with_no_index(site_repo, tmp_path):
    """A snapshot missing `index.html` is refused rather than force-pushed.

    The force-push would otherwise replace a live site with nothing. Staging is
    sabotaged here the only way that can produce an empty snapshot -- dropping the
    `--force` that overrides the exclude file -- to prove the guard is reachable and
    that the published branch survives.

    """
    proj, remote = site_repo
    assert run_script(proj).returncode == 0
    published_before = published_files(remote)

    sabotaged = tmp_path / "sabotaged.sh"
    sabotaged.write_text(
        SCRIPT.read_text().replace("add --all --force .", "add --all .")
    )
    (proj / ".git" / "info" / "exclude").write_text("*\n")

    result = run_script(proj, script=sabotaged)
    assert result.returncode != 0
    assert "no top-level index.html" in result.stderr
    assert published_files(remote) == published_before


def test_refuses_the_remotes_default_branch(site_repo):
    """Publishing over the remote's default branch is refused."""
    proj, _ = site_repo
    result = run_script(proj, {"PUBLISH_DOCS_GH_PAGES_BRANCH": "main"})
    assert result.returncode != 0
    assert "default branch 'main'" in result.stderr


@pytest.mark.parametrize("branch", ["master", "develop"])
def test_refuses_reserved_branch_names(site_repo, branch):
    """Publishing to a branch that is conventionally a source branch is refused."""
    proj, _ = site_repo
    result = run_script(proj, {"PUBLISH_DOCS_GH_PAGES_BRANCH": branch})
    assert result.returncode != 0
    assert f"Refusing to publish to '{branch}'" in result.stderr


def test_requires_an_index_html_in_the_site_dir(site_repo):
    """A site directory with no top-level `index.html` is rejected up front."""
    proj, _ = site_repo
    (proj / "results" / "docs" / "index.html").unlink()
    result = run_script(proj)
    assert result.returncode != 0
    assert "Missing index.html" in result.stderr


def test_requires_the_site_dir_to_exist(site_repo):
    """A missing site directory is rejected up front."""
    proj, _ = site_repo
    result = run_script(proj, {"PUBLISH_DOCS_GH_PAGES_SITE_DIR": "results/nosuch"})
    assert result.returncode != 0
    assert "Directory not found" in result.stderr


def test_requires_git_identity(site_repo):
    """Without `user.name`/`user.email` the script stops before touching the remote."""
    proj, _ = site_repo
    git("config", "--unset", "user.email", cwd=proj)
    result = run_script(proj)
    assert result.returncode != 0
    assert "user.email is not configured" in result.stderr
