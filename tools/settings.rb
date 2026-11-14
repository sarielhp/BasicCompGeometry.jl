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
  # JET runtime-dispatch check of the hot entry points, in its own environment
  # (test/jet) so Pkg.test, CI and package users never install JET.
  FULL_CHECKS = [
    ['julia', '--project=test/jet', '--startup-file=no', '-e', 'using Pkg; Pkg.instantiate()'],
    ['julia', '--project=test/jet', '--startup-file=no', 'test/jet/jet_opt.jl']
  ].freeze
  VERIFICATION_INPUTS = ['test/Project.toml', 'docs/Project.toml', 'test/jet/Project.toml'].freeze
end
