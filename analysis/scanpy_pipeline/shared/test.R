if (file.exists("/.dockerenv")) {
   docker = TRUE
}else {
   docker = FALSE
}

print(docker)
