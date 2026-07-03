# Releasing DIRT

This document describes how a tagged DIRT release is cut and how each release is
archived to [Zenodo](https://zenodo.org/) to obtain a citeable DOI. A DOI is a
[JOSS](https://joss.theoj.org/) requirement at acceptance, so every release that
we intend to cite must be archived.

DIRT follows [Semantic Versioning](https://semver.org/) and keeps a
[Keep a Changelog](https://keepachangelog.com/)-style [`CHANGELOG.md`](CHANGELOG.md).
While the project is in the `0.1.x` pre-1.0 series, minor internal APIs may change
between releases without a major-version bump.

## Remotes

DIRT lives on two remotes:

- **Gitea** (`192.168.0.170`) is the **primary** development remote. All work
  lands on `main` through reviewed pull requests, and the annotated `vX.Y.Z`
  release tag is cut here.
- **GitHub** (`github.com/SueHeir/dirt`) is a **published downstream mirror** of
  the regression-gated Gitea `main`. It is the remote that Zenodo watches, so the
  public GitHub Release is what triggers DOI minting.

## Cutting a release

1. **Land all release content on Gitea `main`.** Every change is a reviewed PR;
   nothing is pushed to `main` directly.

2. **Update the release metadata** (each via its own reviewed PR):
   - Add a dated `## [X.Y.Z] - YYYY-MM-DD` section to `CHANGELOG.md` describing
     `Added` / `Changed` / `Removed` / `Fixed`, and add the matching
     `[X.Y.Z]: .../compare/vPREV...vX.Y.Z` link at the bottom.
   - Bump `version:` and `date-released:` in [`CITATION.cff`](CITATION.cff) to
     match.

3. **Cut the annotated tag** on the `main` commit whose `CHANGELOG.md` and
   `CITATION.cff` describe `X.Y.Z`:

   ```bash
   git fetch origin
   git tag -a vX.Y.Z <release-commit-on-main> \
     -m "DIRT vX.Y.Z" -m "<one-line summary of the CHANGELOG section>"
   git push origin vX.Y.Z
   ```

   The tag is **annotated** (not lightweight) so it carries an author, date, and
   message. Tags are immutable once published — never move or force-push a
   release tag; cut a new patch version instead.

4. **Publish to GitHub.** Push the regression-gated Gitea `main` (and a dated
   stable tag) to the public GitHub mirror:

   ```bash
   ~/projects/automation/bin/publish-to-github.sh
   ```

   This rewrites inter-crate dependency URLs from Gitea to `github.com` and
   pushes `main` for `grass`, `soil`, and `dirt` in dependency order. Push the
   `vX.Y.Z` tag to GitHub as well so it is available as a GitHub Release target.

## GitHub → Zenodo archival (DOI)

Zenodo mints a DOI automatically when a **new GitHub Release** is published on a
repository it is watching. The wiring is:

1. **Enable the integration (maintainer, one-time manual step — not
   automatable).** A maintainer with admin rights on the GitHub repository logs
   in to Zenodo with GitHub, authorizes the Zenodo GitHub OAuth app, and flips
   the toggle **on** for `SueHeir/dirt` at
   <https://zenodo.org/account/settings/github/>. This grants Zenodo a webhook on
   the repo. **This requires an interactive human OAuth grant and cannot be done
   by the fleet or any script.**

2. **Publish a GitHub Release** for the pushed `vX.Y.Z` tag (GitHub → *Releases*
   → *Draft a new release* → choose tag `vX.Y.Z`). When the release is published,
   the Zenodo webhook fires, Zenodo ingests the tarball, and it mints:
   - a **version DOI** unique to `vX.Y.Z`, and
   - a **concept DOI** that always resolves to the latest version.

3. **Record the DOI.** Once Zenodo reports the DOI:
   - Add the Zenodo DOI badge to [`README.md`](README.md).
   - Add a `doi:` (or `identifiers:`) entry to [`CITATION.cff`](CITATION.cff).
   - Reference the DOI in the JOSS submission.

> **Status:** As of this document, the GitHub↔Zenodo integration has **not** been
> enabled and **no DOI has been minted**. Do not cite a DOI until a maintainer
> completes step 1 above and a GitHub Release produces one.

## Checklist

- [ ] `CHANGELOG.md` has a dated `[X.Y.Z]` section + compare link
- [ ] `CITATION.cff` `version` / `date-released` bumped to match
- [ ] Annotated tag `vX.Y.Z` pushed to Gitea `origin`
- [ ] Gitea `main` published to the GitHub mirror; `vX.Y.Z` tag pushed to GitHub
- [ ] (maintainer) GitHub↔Zenodo integration enabled once, before the release
- [ ] GitHub Release published for `vX.Y.Z` → Zenodo DOI minted
- [ ] DOI recorded in `README.md` + `CITATION.cff` + JOSS submission
