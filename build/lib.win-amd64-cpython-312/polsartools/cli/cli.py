# polsartools/cli.py

import click
from .polsar.fp import grvi, rvifp  # Import other functions as needed

@click.group()
def cli():
    """A command line interface for polsartools."""
    pass

@cli.command()
@click.argument('input')  # Add your arguments
@click.argument('output')
def grvi_command(input, output):
    """Run the grvi function."""
    result = grvi(input)  # Call the grvi function
    # Save or print the result as needed
    click.echo(f"Result: {result}")

@cli.command()
@click.argument('input')
@click.argument('output')
def rvifp_command(input, output):
    """Run the rvifp function."""
    result = rvifp(input)
    click.echo(f"Result: {result}")

if __name__ == '__main__':
    cli()
