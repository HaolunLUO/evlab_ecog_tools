function dirFiles = cleanFileNames(dirFiles)
dirFiles = dirFiles(~startsWith({dirFiles.name}, '.'));
end