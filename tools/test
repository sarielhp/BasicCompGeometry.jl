#!/usr/bin/env ruby
# frozen_string_literal: true

require_relative 'tooling'

HELP = <<~TEXT
  Usage: tools/test

  Run the package's test suite: Pkg.test() in the project environment.

  Options:
    -h, --help  Show this help; nothing is run.
TEXT

Tooling.entry(HELP) do
  abort 'Usage: tools/test' unless ARGV.empty?
  ProjectTmp.with_test_env { Language.test }
end
