rule hello_world:
	output:
		"Hello/World.txt",
	shell:
		"""
		echo "Hello, World!" > {output}
		"""
