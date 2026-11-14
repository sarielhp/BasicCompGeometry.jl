# frozen_string_literal: true

module ToolSettings
  # Name of the project's temporary-file tree, <tmpdir>/<user>/<name>/; nil uses the
  # repository directory's name.
  PROJECT_NAME = nil
  TOOLS = {}.freeze
  ARTIFACTS = [].freeze # A package has none; list executable scripts (e.g. 'bin/tool') to check and install them.
  INSTALL_TARGETS = {}.freeze # Artifact => absolute, authorized destination.
  VERSION_FILE = 'Project.toml'
  VERSION_UPDATE = [].freeze
  # A package registered in General is tagged by TagBot once tools/publish's
  # registration merges, so tools/bump pushes no tag. Set true for a project
  # that is never registered.
  RELEASE_TAG = false
  GITHUB_REPO = nil # owner/Name.jl for tools/publish; nil reads it from the upstream remote.
  PUBLISH_WORKFLOW = 'CI.yml' # The workflow whose run on the commit must pass before registering.
  FAST_CHECKS = [].freeze
  # Slower checks, run by tools/gate --full (tools/release, review-cycle):
  # - property/reference tests on seeded random inputs (test/properties);
  # - allocation baselines of the hot entry points (test/alloc);
  # - JET runtime-dispatch check, in its own environment (test/jet) so
  #   Pkg.test, CI and package users never install JET;
  # - the Documenter build (docs/).
  FULL_CHECKS = [
    ['julia', '--project=.', '--startup-file=no', '-e', 'using Test; include("test/properties/properties.jl")'],
    ['julia', '--project=.', '--startup-file=no', 'test/alloc/alloc_check.jl'],
    ['julia', '--project=test/jet', '--startup-file=no', '-e', 'using Pkg; Pkg.instantiate()'],
    ['julia', '--project=test/jet', '--startup-file=no', 'test/jet/jet_opt.jl'],
    # Documenter build: catches broken @ref links, failing @example blocks, and
    # docstrings detached from their function.
    ['julia', '--project=docs', '--startup-file=no', '-e', 'using Pkg; Pkg.instantiate()'],
    ['julia', '--project=docs', '--startup-file=no', 'docs/make.jl']
  ].freeze
  VERIFICATION_INPUTS = ['test/Project.toml', 'docs/Project.toml', 'test/jet/Project.toml'].freeze
end
