package edu.duke.cs.osprey.structure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;

/**
 * This class implements a collection of Residues with fast lookups.
 */
/*
	Believe it or not, residue lookups by string actually became a performance bottleneck.
	
	Sadly, to maintain compatibility with Molecule.residues, we have to subclass ArrayList<Residue>.
	And since the JVM ArrayList implementation has optimized just about every mutator, we have to override
	all those methods here to keep the hash map consistent with the list state.
	
	Oi... what a pain in the butt...
*/
public class Residues extends ArrayList<Residue> {
	
	private static final long serialVersionUID = 4950265664624559802L;
	
	private Map<String,Integer> indicesByNum;
	
	public Residues() {
		indicesByNum = new HashMap<>();
	}
	
	public Residues(Residue ... residues) {
		this();
		addAll(Arrays.asList(residues));
	}
	
	@Override
	public boolean add(Residue res) {
		int index = size();
		super.add(null);
		set(index, res);
		return true;
	}
	
	@Override
	public void add(int index, Residue res) {
		super.add(index, null);
		set(index, res);
	}
	
	@Override
	public Residue set(int index, Residue res) {
		String num = res.getPDBResNumber();
		if (indicesByNum.containsKey(num)) {
			throw new IllegalArgumentException("residue number " + num + " is already present in this collection");
		}
		indicesByNum.put(num, index);
		return super.set(index, res);
	}
	
	@Override
	public boolean remove(Object obj) {
		return remove((Residue)obj);
	}
	
	public boolean remove(Residue res) {
		return remove(res.getPDBResNumber()) != null;
	}
	
	public Residue remove(String num) {
		Integer index = findIndex(num);
		if (index != null) {
			return remove((int)index);
		}
		return null;
	}
	
	@Override
	public Residue remove(int index) {
		Residue res = super.remove(index);
		if (res != null) {
			reindex();
		}
		return res;
	}
	
	@Override
	public void clear() {
		super.clear();
		reindex();
	}
	
	private void reindex() {
		indicesByNum.clear();
		for (int i=0; i<size(); i++) {
			indicesByNum.put(get(i).getPDBResNumber(), i);
		}
	}
	
	@Override
	public boolean contains(Object obj) {
		return super.contains((Residue)obj);
	}

	public boolean contains(String resNum) {
		return getOrNull(resNum) != null;
	}

	@Override
	public int indexOf(Object obj) {
		return indexOf((Residue)obj);
	}
	
	public int indexOf(Residue res) {
		Integer index = findIndex(res);
		if (index != null) {
			return index;
		}
		return -1;
	}
	
	@Override
	public int lastIndexOf(Object obj) {
		return lastIndexOf((Residue)obj);
	}
	
	public int lastIndexOf(Residue res) {
		// residues are unique, so only one exists in each collection
		// forward index is the same as reverse index
		return indexOf(res);
	}
	
	public Integer findIndex(Residue res) {
		return findIndex(res.getPDBResNumber());
	}
	
	public Integer findIndex(String num) {
		return indicesByNum.get(num);
	}
	
	public int findIndexOrThrow(Residue res) {
		return findIndexOrThrow(res.getPDBResNumber());
	}
	
	public int findIndexOrThrow(String num) {
		Integer index = findIndex(num);
		if (index != null) {
			return index;
		}
		throw new NoSuchElementException("no residue with number " + num + " was found");
	}
	
	public Residue getOrThrow(String num) {
		return get(findIndexOrThrow(num));
	}
	
	public Residue getOrNull(String num) {
		Integer index = findIndex(num);
		if (index != null) {
			return get(index);
		}
		return null;
	}
	
	@Override
	public Object clone() {
		Residues residues = (Residues)super.clone();
		residues.indicesByNum = new HashMap<>(this.indicesByNum);
		return residues;
	}
	
	@Override
	public boolean addAll(Collection<? extends Residue> c) {
		return reindexIfChanged(super.addAll(c));
	}
	
	@Override
	public boolean addAll(int index, Collection<? extends Residue> c) {
		return reindexIfChanged(super.addAll(index, c));
	}
	
	@Override
	protected void removeRange(int fromIndex, int toIndex) {
		super.removeRange(fromIndex, toIndex);
		reindex();
	}
	
	@Override
	public boolean removeAll(Collection<?> c) {
		return reindexIfChanged(super.removeAll(c));
	}
	
	@Override
	public boolean retainAll(Collection<?> c) {
		return reindexIfChanged(super.retainAll(c));
	}
	
	@Override
	public boolean removeIf(Predicate<? super Residue> filter) {
		return reindexIfChanged(super.removeIf(filter));
	}
	
	@Override
	public void replaceAll(UnaryOperator<Residue> operator) {
		super.replaceAll(operator);
		reindex();
	}
	
	@Override
	public void sort(Comparator<? super Residue> c) {
		super.sort(c);
		reindex();
	}
	
	private boolean reindexIfChanged(boolean changed) {
		if (changed) {
			reindex();
		}
		return changed;
	}
}
