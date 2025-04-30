from cadetrdm import process_example, Options

if __name__ == "__main__":
    options = Options()
    options.commit_message = "Run large optimizations"
    options.debug = False
    options.push = False
    options.source_directory = "src"
    process_example(options)
