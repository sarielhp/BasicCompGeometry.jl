#!/usr/bin/env ruby
# Snapshot: add all files, commit with timestamp, and push (no version bump).

repo = File.join(__dir__, '..')
timestamp = Time.now.strftime('%Y-%m-%d %H:%M:%S')

system('git', '-C', repo, 'add', '-A') || abort('git add failed')
system('git', '-C', repo, 'commit', '-m', "Snapshot #{timestamp}") || abort('git commit failed')
system('git', '-C', repo, 'push') || abort('git push failed')

puts "Snapshot committed and pushed at #{timestamp}"
