#!/usr/bin/env ruby
# frozen_string_literal: true

require_relative 'tooling'

HELP = <<~TEXT
  Usage: tools/gate [--full]

  Check the tree without changing it: tools/lint, FAST_CHECKS and tools/test.

  Options:
    --full      The release gate: adds the deep review and FULL_CHECKS.
    -h, --help  Show this help; nothing is run.
TEXT

Tooling.entry(HELP) do
  abort 'Usage: tools/gate [--full]' unless ARGV.empty? || ARGV == ['--full']
  Tooling.gate(ARGV == ['--full'])
end
