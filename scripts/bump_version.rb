#!/usr/bin/env ruby
# Bump patch version (0.0.1), commit, push, and tag on GitHub.

PROJECT_FILE = File.join(__dir__, '..', 'Project.toml')

content = File.read(PROJECT_FILE)
old_ver = content[/^version\s*=\s*"(.+?)"/, 1]
abort 'Could not find version in Project.toml' unless old_ver

major, minor, patch = old_ver.split('.').map(&:to_i)
new_ver = "#{major}.#{minor}.#{patch + 1}"

content.sub!(/^version\s*=\s*".*?"/, "version = \"#{new_ver}\"")
File.write(PROJECT_FILE, content)

puts "Bumped version: #{old_ver} -> #{new_ver}"

repo = File.join(__dir__, '..')
system('git', '-C', repo, 'add', 'Project.toml') || abort('git add failed')
system('git', '-C', repo, 'commit', '-m', "Bump version to #{new_ver}") || abort('git commit failed')
system('git', '-C', repo, 'push') || abort('git push failed')
system('git', '-C', repo, 'tag', '-a', "v#{new_ver}", '-m', "Version #{new_ver}") || abort('git tag failed')
system('git', '-C', repo, 'push', '--tags') || abort('git push --tags failed')

puts "Tagged and pushed v#{new_ver}"
