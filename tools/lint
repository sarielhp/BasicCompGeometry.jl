#!/usr/bin/env ruby
# frozen_string_literal: true

require_relative 'tooling'

HELP = <<~TEXT
  Usage: tools/lint

  Static feedback only: julia-audit --strict and julia-static-analysis
  --syntax --format. No tests.

  Options:
    -h, --help  Show this help; nothing is run.
TEXT

Tooling.entry(HELP) do
  abort 'Usage: tools/lint' unless ARGV.empty?
  Language.lint
end
