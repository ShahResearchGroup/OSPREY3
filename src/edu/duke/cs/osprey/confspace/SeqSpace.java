package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.tools.MathTools;

import java.util.*;
import java.util.stream.Collectors;


/**
 * analogous to the conf space, but just for sequences
 */
public class SeqSpace {

	/**
	 * Return the unique sequence space common to all conformation spaces,
	 * or throw an exception.
	 */
	public static SeqSpace reduce(List<SimpleConfSpace> confSpaces) {
		return confSpaces.stream()
			.map(confSpace -> confSpace.seqSpace)
			.reduce((a, b) -> {
				if (!a.equals(b)) {
					throw new IllegalArgumentException(String.format(
						"conf spaces have different sequence spaces:\n%s\n%s",
						a, b
					));
				}
				return a;
			})
			.orElseThrow(() -> new IllegalArgumentException("conf spaces list is empty"));
	}

	/**
	 * Return the sequence space from the conformation space with the greatest number of mutable positions
	 */
	public static SeqSpace max(List<SimpleConfSpace> confSpaces) {
		return confSpaces.stream()
			.map(confSpace -> confSpace.seqSpace)
			.max(Comparator.comparing(seqSpace -> seqSpace.positions.size()))
			.orElseThrow(() -> new IllegalArgumentException("conf spaces list is empty"));
	}

	public class Position implements Comparable<Position> {

		public final SeqSpace seqSpace = SeqSpace.this;

		public final int index;
		public final String resNum;
		public final ResType wildType;
		public final List<ResType> resTypes;
		public final List<ResType> mutations;

		private final Map<String,ResType> resTypesByName;

		private Position(int index, String resNum, String wildType, List<String> resTypes) {

			this.index = index;
			this.resNum = resNum;

			// make the res types
			this.resTypes = new ArrayList<>(resTypes.size());
			for (String resType : resTypes) {
				this.resTypes.add(new ResType(
					this,
					this.resTypes.size(),
					resType
				));
			}

			// index them by name
			resTypesByName = new HashMap<>();
			for (ResType rt : this.resTypes) {
				resTypesByName.put(rt.name, rt);
			}

			// get the wild type
			this.wildType = getResType(wildType);

			// gather the mutants
			mutations = this.resTypes.stream()
				.filter(rt -> rt != this.wildType)
				.collect(Collectors.toList());
		}

		public ResType getResType(String name) {
			return resTypesByName.get(name);
		}

		public ResType getResTypeOrThrow(String name) {
			ResType rt = getResType(name);
			if (rt != null) {
				return rt;
			}
			throw new NoSuchElementException("Res type " + name + " not allowed at position " + resNum + ". Try one of " + resTypes);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof Position && equals((Position)other);
		}

		public boolean equals(Position other) {
			return this.index == other.index
				&& this.resNum.equals(other.resNum)
				&& this.wildType.equals(other.wildType)
				&& this.resTypes.equals(other.resTypes);
		}

		@Override
		public String toString() {
			return String.format("%d:%s", index, resNum);
		}

		@Override
		public int compareTo(Position other) {
			return this.index - other.index;
		}
	}

	public class ResType implements Comparable<ResType> {

		public final Position pos;
		public final int index;
		public final String name;

		public ResType(Position pos, int index, String name) {
			this.pos = pos;
			this.index = index;
			this.name = name;
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof ResType && equals((ResType)other);
		}

		public boolean equals(ResType other) {
			return this.index == other.index
				&& this.name.equals(other.name);
		}

		@Override
		public String toString() {
			return String.format("%d:%s", index, name);
		}

		@Override
		public int compareTo(ResType other) {
			return this.index - other.index;
		}

		public boolean isWildType() {
			return pos.wildType == this;
		}

		public boolean isMutation() {
			return pos.wildType != this;
		}

		/**
		 * Returns the residue type name in lower case if wild-type, upper case if mutation.
		 */
		public String mutationName() {
			if (isWildType()) {
				return name.toLowerCase();
			} else {
				return name.toUpperCase();
			}
		}
	}

	public final List<Position> positions;

	private final Map<String,Position> positionsByResNum;

	public SeqSpace(SimpleConfSpace confSpace) {

		// make the positions
		positions = new ArrayList<>(confSpace.mutablePositions.size());
		for (SimpleConfSpace.Position pos : confSpace.mutablePositions) {
			positions.add(new Position(
				positions.size(),
				pos.resNum,
				pos.resFlex.wildType,
				pos.resTypes
			));
		}

		// index by res num
		positionsByResNum = new HashMap<>();
		for (Position pos : positions) {
			positionsByResNum.put(pos.resNum, pos);
		}
	}

	public List<String> getResNums() {
		return positions.stream()
			.map(pos -> pos.resNum)
			.collect(Collectors.toList());
	}

	public Position getPosition(String resNum) {
		return positionsByResNum.get(resNum);
	}

	public Position getPositionOrThrow(String resNum) {
		Position pos = getPosition(resNum);
		if (pos != null) {
			return pos;
		}
		throw new NoSuchElementException("no pos found in this sequence space at residue " + resNum + ". try one of " + getResNums());
	}

	public Sequence makeUnassignedSequence() {
		return new Sequence(this);
	}

	public Sequence makeWildTypeSequence() {
		Sequence seq = new Sequence(this);
		seq.fillWildType();
		return seq;
	}

	/**
	 * Make a sequence with the given residue type names
	 *
	 * residue type names must be given in the same order as the positions in the sequence space
	 */
	public Sequence makeSequence(List<String> resTypes) {

		// just in case...
		if (resTypes.size() != positions.size()) {
			throw new IllegalArgumentException(String.format("expected %d residue types, but only got %d: %s",
				positions.size(),
				resTypes.size(),
				resTypes
			));
		}

		Sequence seq = makeUnassignedSequence();
		for (int i=0; i<positions.size(); i++) {
			Position pos = positions.get(i);
			seq.set(pos, resTypes.get(i));
		}
		return seq;
	}

	public List<Sequence> getMutants() {
		return getMutants(positions.size());
	}

	public List<Sequence> getMutants(int maxSimultaneousMutations) {
		return getMutants(maxSimultaneousMutations, false);
	}

	public List<Sequence> getMutants(int maxSimultaneousMutations, boolean reversePositionOrder) {

		List<Sequence> sequences = new ArrayList<>();

		// get all possible combinations of mutations
		List<List<SeqSpace.Position>> powersetOfPositions = MathTools.powersetUpTo(positions, maxSimultaneousMutations);

		// reverse to match order of old K* code if needed
		if (reversePositionOrder) {
			Collections.reverse(powersetOfPositions);
		}

		for (List<SeqSpace.Position> positions : powersetOfPositions) {

			// collect the mutations (res types except for wild type) for these positions into a simple list list
			List<List<ResType>> mutationsByPos = positions.stream()
				.map(pos -> pos.mutations)
				.collect(Collectors.toList());

			// enumerate all the combinations of res types
			for (List<ResType> mutations : MathTools.cartesianProduct(mutationsByPos)) {

				// build the complex sequence
				Sequence sequence = makeWildTypeSequence();
				for (ResType rt : mutations) {
					sequence.set(rt.pos, rt);
				}
				sequences.add(sequence);
			}
		}

		return sequences;
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append("Residue Types:");
		for (SeqSpace.Position pos : positions) {
			buf.append("\n\t");
			buf.append(pos.index);
			buf.append(":");
			buf.append(pos.resNum);
			buf.append("  [");
			for (ResType rt : pos.resTypes) {
				buf.append(" ");
				buf.append(rt.index);
				buf.append(":");
				if (pos.wildType == rt) {
					buf.append(rt.name.toLowerCase());
				} else {
					buf.append(rt.name);
				}
			}
			buf.append(" ]");
		}
		return buf.toString();
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof SeqSpace && equals((SeqSpace)other);
	}

	public boolean equals(SeqSpace other) {
		return this.positions.equals(other.positions);
	}
}