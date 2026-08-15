#!/usr/bin/env ruby
# frozen_string_literal: true

require 'zlib'
require 'optparse'

#
# verify_page_widths.rb
#
# Extracts content stream of a given page from output/prophet.pdf,
# decompresses the FlateDecode stream, parses the PDF graphics state operators
# (including CTM matrix `cm`, line width `w`, path constructions `m`/`l`, and strokes `S`/`s`),
# calculates exact device point widths for all drawn segments, and validates the statistical distribution.
#

def extract_pdf_page_streams(pdf_path)
  content = File.binread(pdf_path)
  streams = []

  content.scan(/(\d+)\s+(\d+)\s+obj\s*<<(.*?)>>\s*stream\r?\n(.*?)\r?\nendstream/m) do |id, _gen, dict, raw|
    if dict.include?('/FlateDecode')
      begin
        decomp = Zlib::Inflate.inflate(raw)
        streams << { id: id.to_i, dict: dict, data: decomp }
      rescue Zlib::DataError, Zlib::BufError
        # Try raw deflate without zlib header
        begin
          z = Zlib::Inflate.new(-Zlib::MAX_WBITS)
          decomp = z.inflate(raw)
          z.close
          streams << { id: id.to_i, dict: dict, data: decomp }
        rescue StandardError
          # Ignore non-decodable streams
        end
      end
    else
      streams << { id: id.to_i, dict: dict, data: raw }
    end
  end

  # Filter and sort page content streams by page title
  page_streams = []
  streams.each do |s|
    if (m = s[:data].match(/Page\s+(\d+)\s*\/\s*(\d+)/))
      page_num = m[1].to_i
      page_streams << { page: page_num, id: s[:id], data: s[:data], title: m[0] }
    elsif s[:data].include?('Verification')
      # Extra verification pages
      page_streams << { page: nil, id: s[:id], data: s[:data], title: 'Verification Page' }
    end
  end

  page_streams
end

def analyze_page_stroke_widths(stream_data)
  # Look for CTM scaling matrix (e.g. `700 0 0 -700 50 750 cm`)
  ctm_scale = 1.0
  if (cm_match = stream_data.match(/([0-9\.\-]+)\s+0\s+0\s+([0-9\.\-]+)\s+[0-9\.\-]+\s+[0-9\.\-]+\s+cm/))
    sx = cm_match[1].to_f.abs
    sy = cm_match[2].to_f.abs
    ctm_scale = (sx + sy) / 2.0 if sx.positive?
  end

  # Parse graphics operators in sequence
  tokens = stream_data.scan(/[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?|[a-zA-Z*]+|\[|\]/)

  current_width_user = nil
  stroke_widths_pt = []
  segments = []
  current_path_points = []

  i = 0
  while i < tokens.length
    tok = tokens[i]

    if tok == 'w' && i.positive?
      current_width_user = tokens[i - 1].to_f
    elsif tok == 'm' && i >= 2
      # Start new subpath
      x = tokens[i - 2].to_f
      y = tokens[i - 1].to_f
      current_path_points = [[x, y]]
    elsif tok == 'l' && i >= 2
      # Line to
      x = tokens[i - 2].to_f
      y = tokens[i - 1].to_f
      if current_path_points.any?
        p_prev = current_path_points.last
        effective_width = current_width_user ? (current_width_user * ctm_scale) : 1.0
        segments << { from: p_prev, to: [x, y], width_pt: effective_width }
      end
      current_path_points << [x, y]
    elsif %w[S s B b B* b*].include?(tok)
      if current_width_user
        effective_width = current_width_user * ctm_scale
        stroke_widths_pt << effective_width
      end
      current_path_points = []
    elsif tok == 'n' || tok == 'f' || tok == 'F' || tok == 'f*'
      current_path_points = []
    end

    i += 1
  end

  {
    ctm_scale: ctm_scale,
    stroke_widths: stroke_widths_pt,
    segments: segments
  }
end

def main
  target_page_num = 100
  pdf_path = File.expand_path('../output/prophet.pdf', __dir__)

  OptionParser.new do |opts|
    opts.banner = "Usage: ruby scripts/verify_page_widths.rb [options]"
    opts.on("-p", "--page NUM", Integer, "Page number to analyze (default: #{target_page_num})") do |p|
      target_page_num = p
    end
    opts.on("-f", "--file PATH", String, "Path to PDF file (default: #{pdf_path})") do |f|
      pdf_path = f
    end
  end.parse!

  unless File.exist?(pdf_path)
    puts "Error: PDF file not found at #{pdf_path}"
    exit 1
  end

  page_streams = extract_pdf_page_streams(pdf_path)

  puts "=" * 76
  puts "PDF Segment Stroke Width Verification Tool"
  puts "=" * 76
  puts "PDF File:        #{pdf_path}"
  puts "Identified Pages: #{page_streams.length}"
  puts "Target Page:     #{target_page_num}"

  target_stream = page_streams.find { |s| s[:page] == target_page_num }

  unless target_stream
    puts "Error: Could not locate content stream for Page #{target_page_num}."
    puts "Available numbered pages: #{page_streams.map { |s| s[:page] }.compact.sort.inspect}"
    exit 1
  end

  analysis = analyze_page_stroke_widths(target_stream[:data])
  widths = analysis[:stroke_widths]
  segments = analysis[:segments]
  ctm_scale = analysis[:ctm_scale]

  puts "\nStream Analysis for Page #{target_page_num}:"
  puts "  Decompressed stream size: #{target_stream[:data].bytesize} bytes"
  puts "  Detected CTM scale factor: #{ctm_scale} (converts user units -> device points)"
  puts "  Total stroked paths:      #{widths.length}"
  puts "  Total line segments:       #{segments.length}"

  if widths.empty?
    puts "Error: No strokes found on Page #{target_page_num}."
    exit 1
  end

  # Statistical computations
  min_w = widths.min
  max_w = widths.max
  mean_w = widths.sum / widths.length.to_f
  sorted = widths.sort
  median_w = sorted[sorted.length / 2]
  variance = widths.map { |w| (w - mean_w)**2 }.sum / widths.length
  std_dev = Math.sqrt(variance)
  unique_widths = widths.map { |w| w.round(1) }.uniq.sort

  # Theoretical values for Uniform(1, 128):
  # Mean = (1 + 128)/2 = 64.5
  # StdDev = sqrt((128 - 1)^2 / 12) = sqrt(16129 / 12) = 36.66
  theoretical_mean = 64.5
  theoretical_std = 36.66

  puts "\n" + "-" * 76
  puts "Statistical Verification of Widths on Page #{target_page_num} (Expected: Uniform in [1, 128]):"
  puts "-" * 76
  printf("  Strokes Sample Count: %d\n", widths.length)
  printf("  Unique Width Values:  %d distinct width levels\n", unique_widths.length)
  printf("  Minimum Stroke Width: %.2f pt  (Expected min: ~1 pt)\n", min_w)
  printf("  Maximum Stroke Width: %.2f pt  (Expected max: ~128 pt)\n", max_w)
  printf("  Sample Mean Width:    %.2f pt  (Theoretical: %.2f pt, Error: %.2f%%)\n",
         mean_w, theoretical_mean, ((mean_w - theoretical_mean).abs / theoretical_mean) * 100)
  printf("  Sample Median Width:  %.2f pt  (Theoretical: %.2f pt)\n", median_w, theoretical_mean)
  printf("  Sample Std Deviation: %.2f pt  (Theoretical: %.2f pt, Error: %.2f%%)\n",
         std_dev, theoretical_std, ((std_dev - theoretical_std).abs / theoretical_std) * 100)

  puts "\nDistribution Histogram across 8 Equal Bins in [1, 128] pt:"
  puts "-" * 76
  bins = [
    (1..16),
    (17..32),
    (33..48),
    (49..64),
    (65..80),
    (81..96),
    (97..112),
    (113..128)
  ]

  bins.each do |range|
    count = widths.count { |w| range.cover?(w.round) }
    pct = (count.to_f / widths.length) * 100
    bar = '#' * (pct * 0.8).round
    printf("  Bin [%3d - %3d pt]: %4d segments (%5.1f%%) | %s\n", range.first, range.last, count, pct, bar)
  end

  puts "\nFirst 15 individual segment strokes drawn on Page #{target_page_num}:"
  segments.first(15).each_with_index do |seg, idx|
    p1 = seg[:from].map { |c| sprintf('%.3f', c) }
    p2 = seg[:to].map { |c| sprintf('%.3f', c) }
    printf("  Segment %2d: (%s, %s) -> (%s, %s) | Width = %6.2f pt\n",
           idx + 1, p1[0], p1[1], p2[0], p2[1], seg[:width_pt])
  end

  puts "\n" + "=" * 76
  # Criteria for pass:
  # 1. Broad support covering close to 1 and 128
  # 2. Many unique values (> 50)
  # 3. Mean within 15% of 64.5
  pass = unique_widths.length >= 50 && min_w <= 10.0 && max_w >= 115.0 && (mean_w - theoretical_mean).abs < 10.0
  if pass
    puts "RESULT: [PASSED] - Page #{target_page_num} segments have genuine random widths spanning [1, 128] pt!"
  else
    puts "RESULT: [FAILED] - Width distribution does not conform to Uniform(1, 128)."
  end
  puts "=" * 76
end

main if __FILE__ == $PROGRAM_NAME
