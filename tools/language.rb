# frozen_string_literal: true

module Language
  SECTION = /^\s*\[/
  VERSION_ASSIGNMENT = /^\s*version\s*=\s*"(?<value>[^"\n]+)"\s*(?:#.*)?$/
  PACKAGE_KEYS = /^\s*(name|uuid)\s*=/
  # Direct uses of the system temporary directory (see Tooling.check_temporary_files):
  # a parent directory must be given (mktempdir(parent), mktemp(parent), ...).
  TEMP_CALLS = %r{\b(?:mktempdir|mktemp|tempname)\(\s*[;)]|\btempdir\(\s*\)|"/tmp\b} # tmp-exempt: the pattern itself
  COMMENT = /#.*/

  module_function

  def name
    'Julia'
  end

  def setup
    raise 'Project.toml missing' unless File.file?('Project.toml')
  end

  def temp_policy_files
    Dir.glob('{src,test,ext,bin}/**/*.jl') + Dir.glob('*.jl')
  end

  def test_file?(path)
    path.start_with?('test/')
  end

  def tool_paths
    [File.expand_path('~/.juliaup/bin'), File.expand_path('~/bin')]
  end

  def required_inputs
    ['Project.toml', 'Manifest.toml', 'JuliaProject.toml', '.JuliaFormatter.toml']
  end

  def verification_tools
    %w[julia julia-audit julia-static-analysis]
  end

  def environment
    version = Tooling.query(Tooling.executable('julia'), '--startup-file=no', '--version')
    [version, ENV.to_h.select { |key, _| key.start_with?('JULIA') }]
  end

  def preflight(_full)
    verification_tools.each { |tool| Tooling.executable(tool) }
  end

  # Runs Julia code in the project's own environment. --startup-file=no keeps a
  # personal startup.jl (Revise, OhMyREPL) out of every check.
  def julia(code)
    Tooling.run(Tooling.executable('julia'), '--project=.', '--startup-file=no', '-e', code)
  end

  # A Julia package has no build output; executable scripts are listed in
  # ARTIFACTS so tools/build checks them and tools/install can copy them.
  def artifacts
    ToolSettings::ARTIFACTS
  end

  def build
    outputs = artifacts.map { |path| Tooling.contained_path(path) }
    outputs.each do |path|
      raise "build artifact missing or not executable: #{path}" unless File.file?(path) && File.executable?(path)
    end
    julia('using Pkg; Pkg.precompile()')
    puts 'Build passed'
  end

  def lint
    preflight(false)
    Tooling.run(Tooling.executable('julia-audit'), '--strict', '--quiet')
    Tooling.run(Tooling.executable('julia-static-analysis'), '--syntax', '--format', '.')
    puts 'Lint passed'
  end

  def test
    julia('using Pkg; Pkg.test()')
  end

  def gate_start; end

  def gate_finish
    Tooling.script('build')
  end

  # JET and ExplicitImports are slow and need the tools environment, so they run
  # only in the full gate; Aqua runs inside the test suite.
  def deep_review
    Tooling.run(Tooling.executable('julia-static-analysis'), '--jet', '--imports', '.')
  end

  def package?
    File.foreach(ToolSettings::VERSION_FILE).take_while { |line| !line.match?(SECTION) }.count do |line|
      line.match?(PACKAGE_KEYS)
    end == 2
  end

  # The version is the top-level `version` key: the one before the first table.
  def version_match
    path = ToolSettings::VERSION_FILE
    Tooling.contained_path(path)
    raise "project file missing: #{path}" unless File.file?(path)

    offset = 0
    File.read(path).each_line do |line|
      break if line.match?(SECTION)

      match = line.match(VERSION_ASSIGNMENT)
      return [path, match[:value], offset + match.begin(:value), match[:value].length] if match

      offset += line.length
    end
    raise "no top-level version in #{path}; add version = \"0.1.0\" before the first table"
  end

  def version
    version_match.first(2)
  end

  def update_version(updated)
    path, _, start, length = version_match
    contents = File.read(path)
    contents[start, length] = updated
    File.write(path, contents)
    return unless File.file?('CHANGELOG.md')

    contents = File.read('CHANGELOG.md')
    contents.sub!('## [Unreleased]', "## [Unreleased]\n\n## [#{updated}] - #{Date.today}")
    File.write(Tooling.contained_path('CHANGELOG.md'), contents)
  end

  def synchronize_version; end

  def next_milestone
    0
  end

  def prepare_milestone(_count); end

  def milestone_checks(_count); end
end
