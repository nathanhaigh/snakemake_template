rule hello_world:
	output:
		"{cheer}/{world}.txt",
	shell:
		"""
		echo "{wildcards.cheer}, {wildcards.world}!" > {output}
		"""
