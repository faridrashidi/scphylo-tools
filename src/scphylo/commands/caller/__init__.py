"""Caller Module."""

from collections import OrderedDict

import click

from scphylo.commands.caller._1rename import rename
from scphylo.commands.caller._1sra import sra
from scphylo.commands.caller._2fastp import fastp
from scphylo.commands.caller._3bwa import bwa
from scphylo.commands.caller._3star import star
from scphylo.commands.caller._4gatk import gatk
from scphylo.commands.caller._4rsem import rsem
from scphylo.commands.caller._5hapcaller import hapcaller
from scphylo.commands.caller._5mpileup import mpileup
from scphylo.commands.caller._5mutect2 import mutect2
from scphylo.commands.caller._5sequenza import sequenza
from scphylo.commands.caller._5strelka2 import strelka2
from scphylo.commands.caller._5varscan2 import varscan2
from scphylo.commands.caller._6bamreadcount import bamreadcount
from scphylo.commands.caller._6bamsnap import bamsnap
from scphylo.commands.caller._6defuse import defuse
from scphylo.commands.caller._6pyclone import pyclone
from scphylo.commands.caller._6snpeff import snpeff
from scphylo.commands.caller._6varfilter import varfilter
from scphylo.commands.caller._6vartrix import vartrix


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.

    Make sure you initialize the `self.commands` with OrderedDict instance.
    With decorator, use::
        @click.group(cls=NaturalOrderGroup, commands=OrderedDict())
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.

        If the dict is OrderedDict, it will preserve the order commands
        were added.
        """
        return self.commands.keys()


@click.group(
    name="caller",
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings={"max_content_width": 300, "terminal_width": 300},
)
def cli_caller():
    r"""Caller module.

    FastQ: {sample}_1.fastq.gz or {sample}.fastq.gz\n
    BAM: {sample}.transcript.bam for counting\n
    BAM: {sample}.mapped.bam -> {sample}.markdup.bam\n
    """
    return None


cli_caller.add_command(rename)
cli_caller.add_command(sra)
cli_caller.add_command(fastp)
cli_caller.add_command(bwa)
cli_caller.add_command(star)
cli_caller.add_command(gatk)
cli_caller.add_command(hapcaller)
cli_caller.add_command(mutect2)
cli_caller.add_command(strelka2)
cli_caller.add_command(varscan2)
cli_caller.add_command(mpileup)
cli_caller.add_command(sequenza)
cli_caller.add_command(varfilter)
cli_caller.add_command(snpeff)
cli_caller.add_command(rsem)
cli_caller.add_command(bamreadcount)
cli_caller.add_command(bamsnap)
cli_caller.add_command(defuse)
cli_caller.add_command(vartrix)
cli_caller.add_command(pyclone)
