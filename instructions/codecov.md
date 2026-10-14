# Enable Codecov for BasicCompGeometry.jl

The coverage upload currently fails with `Repository not found`. Package tests
and coverage-file generation succeed; Codecov needs to recognize
`sarielhp/BasicCompGeometry.jl` before it can accept the report.

The repository already uploads `lcov.info` through GitHub OIDC authentication.
This setup uses GitHub-issued credentials rather than a `CODECOV_TOKEN` secret.
See [Codecov's OIDC instructions](https://github.com/codecov/codecov-action#using-oidc).

## 1. Give the Codecov GitHub app access to this repository

1. Sign in to GitHub as **sarielhp**.
2. Open the [Codecov GitHub app](https://github.com/apps/codecov).
3. Click **Install** and choose the **sarielhp** personal account.
4. Under repository access, choose **Only select repositories**.
5. Select **BasicCompGeometry.jl** and click **Install**.

If Codecov is already installed, open
[GitHub's Installed GitHub Apps settings](https://github.com/settings/installations),
find **Codecov**, and click **Configure**. Add **BasicCompGeometry.jl** to the
selected repositories and click **Save**. Keep any previously selected
repositories selected.

GitHub describes repository selection in its
[app installation guide](https://docs.github.com/en/apps/using-github-apps/installing-a-github-app-from-a-third-party).

## 2. Connect your account and select the repository in Codecov

1. Open [Codecov for sarielhp](https://app.codecov.io/gh/sarielhp).
2. Sign in with **GitHub**, using the same **sarielhp** account, and complete
   the authorization prompt if one appears.
3. Select **sarielhp** as the account/owner in Codecov.
4. Find **BasicCompGeometry.jl** in the repository list. Check **Not Configured**
   if it does not appear in the configured list.
5. Open the repository and click **Setup repo** if offered. Complete the
   repository onboarding until you can open its dashboard or configuration.
6. Under **Configuration → General**, check that the default branch is **main**.

Codecov may show instructions to copy an upload token. Our workflow already
authenticates with OIDC, so retain the existing workflow configuration. A blank
coverage dashboard is expected until the first upload succeeds.

If the repository is missing, return to the GitHub app configuration and confirm
its access to **BasicCompGeometry.jl**. Allow a few moments for synchronization,
then refresh Codecov. These account and repository steps follow
[Codecov's quick start](https://docs.codecov.com/docs/quick-start).

## 3. Retry the coverage job

1. Open the repository's
   [Documentation workflow runs](https://github.com/sarielhp/BasicCompGeometry.jl/actions/workflows/Documentation.yml).
2. Select the most recent completed run on **main** that contains a separate
   **coverage** job.
3. Open **coverage** and confirm its failed step is **Publish coverage to Codecov**.
4. At the top right of the run page, select **Re-run jobs → Re-run failed jobs**,
   then confirm **Re-run jobs**.

If coverage is the only failed job, this retries the upload using the saved
coverage artifact. GitHub documents these controls in its
[rerun guide](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/re-run-workflows-and-jobs).

Alternatively, use the GitHub CLI. First list runs:

```bash
gh run list --repo sarielhp/BasicCompGeometry.jl \
  --workflow Documentation.yml --branch main --limit 5
```

Replace `RUN_ID` below with the ID of the selected completed run:

```bash
gh run rerun RUN_ID --repo sarielhp/BasicCompGeometry.jl --failed
gh run watch RUN_ID --repo sarielhp/BasicCompGeometry.jl --exit-status
```

You can also tell Codex that the repository is enabled, and it can retry and
verify the upload.

## 4. Verify the result

1. The **coverage** job and **Publish coverage to Codecov** step should turn green.
2. Open the
   [repository's Codecov dashboard](https://app.codecov.io/gh/sarielhp/BasicCompGeometry.jl)
   and confirm that **main** has a coverage percentage and the uploaded commit.
   Processing may take a few moments after the upload.
3. Refresh the [GitHub README](https://github.com/sarielhp/BasicCompGeometry.jl#readme).
   Its coverage badge should display the percentage once the report is processed.

If the upload still says `Repository not found`, check both the GitHub app's
selected repositories and the Codecov account/owner. If it reports a different
error, provide the **Publish coverage to Codecov** log for diagnosis.

Enabling Codecov and rerunning the coverage job do not request Julia registration.
