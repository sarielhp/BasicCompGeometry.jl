#!/usr/bin/env ruby
# frozen_string_literal: true

require 'fileutils'
require 'minitest/autorun'
require 'open3'
require 'rexml/document'
require 'timeout'
require 'yaml'

ROOT = File.expand_path('..', __dir__).freeze
CASES = YAML.safe_load_file(File.join(__dir__, 'example_smoke.yml')).freeze
RUN_DIR = File.join(ROOT, 'output', 'tmp', 'example-smoke', "run-#{Process.pid}").freeze

class ExampleSmokeTest < Minitest::Test
  def test_every_example_has_a_smoke_case
    programs = Dir.glob('**/*.jl', base: File.join(ROOT, 'examples')).reject { |f| f.start_with?('support/') }
    assert_equal programs.sort, CASES.map { |c| c.fetch('program') }.uniq.sort
    assert_equal CASES.size, CASES.map { |c| c.fetch('name') }.uniq.size
    CASES.each { |c| refute_empty c.fetch('purpose') }
  end

  def run_example(spec, directory)
    FileUtils.mkdir_p(directory)
    env = {'QUICKENV_DISABLE_AUTO' => 'true', 'BCG_EXAMPLE_SMOKE' => 'true',
           'BCG_EXAMPLE_OUTPUT_DIR' => directory, 'JULIA_LOAD_PATH' => '@:@stdlib'}
    command = [ENV.fetch('JULIA', 'julia'), '--startup-file=no', "--project=#{ROOT}/examples"]
    if ENV['BCG_EXAMPLE_COVERAGE'] == 'true'
      command.concat(["--code-coverage=@#{ROOT}", "--code-coverage=#{directory}/coverage.info"])
    end
    command.concat([File.join(ROOT, 'examples', spec.fetch('program')), *spec.fetch('args', [])])
    log_path = File.join(directory, 'run.log')
    status = nil
    pid = Process.spawn(env, *command, chdir: ROOT, out: log_path, err: [:child, :out], pgroup: true)
    begin
      status = Timeout.timeout(Integer(ENV.fetch('BCG_EXAMPLE_TIMEOUT', '240'))) { Process.wait2(pid).last }
    ensure
      terminate_example(pid) unless status
    end
    assert status.success?, "#{spec.fetch('name')} failed; log: #{log_path}\n#{File.read(log_path).lines.last(35).join}"
    File.read(log_path)
  end

  def terminate_example(pid)
    Process.kill('KILL', -pid)
    Process.wait(pid)
  rescue Errno::ESRCH, Errno::ECHILD
    # The process may have exited between the timeout and cleanup.
  end

  def check_file(path, expected_pages)
    assert File.file?(path), "Missing output: #{path}"
    assert File.size(path).positive?, "Empty output: #{path}"
    case File.extname(path)
    when '.pdf'
      output, status = Open3.capture2e('qpdf', '--check', path)
      assert status.success?, "Invalid PDF #{path}: #{output}"
      pages, status = Open3.capture2e('qpdf', '--show-npages', path)
      assert status.success?, pages
      assert pages.to_i.positive?, "No pages: #{path}"
      assert_equal expected_pages, pages.to_i, path if expected_pages
    when '.ipe'
      assert_equal 'ipe', REXML::Document.new(File.read(path)).root.name
    when '.png'
      assert_equal "\x89PNG\r\n\x1a\n".b, File.binread(path, 8)
    when '.gif'
      frames, status = Open3.capture2e('ffprobe', '-v', 'error', '-count_frames', '-select_streams', 'v:0',
                                     '-show_entries', 'stream=nb_read_frames', '-of', 'csv=p=0', path)
      assert status.success?, frames
      assert_operator frames.to_i, :>, 1, "Animation has no frames: #{path}"
    when '.html'
      assert_includes File.read(path), '.svg'
      slides = Dir.glob(File.join(File.dirname(path), '**', '*.svg'))
      refute_empty slides
      slides.each { |slide| assert_equal 'svg', REXML::Document.new(File.read(slide)).root.name }
    end
  end

  CASES.each do |spec|
    define_method("test_example_#{spec.fetch('name')}") do
      puts "\n#{spec.fetch('name')}: #{spec.fetch('purpose')}"
      directory = File.join(RUN_DIR, spec.fetch('name'))
      log = run_example(spec, directory)
      spec.fetch('stdout', []).each { |marker| assert_includes log, marker }
      spec.fetch('files', {}).each { |file, pages| check_file(File.join(directory, file), pages) }
    end
  end
end

Minitest.after_run do
  reports = Dir.glob(File.join(RUN_DIR, '*', 'coverage.info')).sort
  next if reports.empty?

  File.open(File.join(RUN_DIR, 'lcov.info'), 'w') do |merged|
    reports.each do |report|
      File.read(report).split("end_of_record\n").each do |record|
        path = record[/^SF:(.+)$/, 1]
        next unless path&.start_with?("#{ROOT}/src/", "#{ROOT}/ext/")

        merged.write(record.sub("SF:#{ROOT}/", 'SF:') + "end_of_record\n")
      end
    end
  end
  puts "Example coverage: #{RUN_DIR}/lcov.info"
end
